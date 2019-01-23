#' Concatenate a date string & suffix
#'
#' @param suffix a character string.
#' @param date_str a character string representing the date.
#' @importFrom glue glue
file_name <- function(suffix, date_str) {
    glue::glue("data/{date_str}_{suffix}.csv")
}

#' Read a csv file with specified date string & suffix
#'
#' @param suffix a character string.
#' @param date_str a character string representing the date.
#' @param ... additional arguments for the column types in readr::read_csv().
#' @importFrom readr read_csv
read_data_csv <- function(suffix, date_str, ...) {
    readr::read_csv(file_name(suffix, date_str), col_types = list(...))
}

#' Extract data from data/<date_str>.json & save as csv files
#'
#' @param date_str a character string representing the date.
#' @import dplyr
#' @importFrom magrittr %>% extract2 is_in subtract %<>%
#' @importFrom purrr map_lgl map_int map_df map map2 map_chr
#' @importFrom lazyeval interp
#' @importFrom rjson fromJSON
#' @importFrom glue glue
#' @importFrom anytime anytime
#' @importFrom tidyr unnest nest
#' @importFrom lubridate floor_date ceiling_date
#' @importFrom tibble data_frame as_tibble
#' @importFrom readr write_csv
#' @export
extraction <- function(date_str) {
    ## 00) prelim
    foo <- function(x) {
        ifelse(purrr::map_lgl(x, is.null), as.character(NA),
        ifelse(purrr::map_int(x, length) == 1L, purrr::map(x, unlist, recursive = F), x))
    }
    bar <- function(df, n) {
        names(which(unlist(dplyr::summarise_all(df, dplyr::funs((function(x) {length(unlist(x, F))})))) != n))
    }
    baz <- function(df) {
        names(which(unlist(dplyr::summarise_all(df, dplyr::funs((function(x) {any(x == "list")}))))))
    }
    qux <- function(vec) {
        ## the survival window function
        1.0 - (rank(vec, "keep") - 0.5) / sum(!is.na(vec))
    }
    quux <- function(l) {
        ## the user non-window function
        df <- dplyr::mutate(tibble::as_tibble(foo(l)), utc_offset = as.numeric(.data$utc_offset), id = as.numeric(.data$id))
        names(df) <- paste("user", names(l), sep = "_")
        df
    }
    ## 01) tweets
    tweets <- rjson::fromJSON(file = glue::glue("data/{date_str}.json")) # raw data, list of lists
    allnames <- unique(unlist(purrr::map(tweets, names))) # all var names
    df0.tweets <- tweets %>% purrr::map_df(purrr::map, list) %>% dplyr::mutate_all(dplyr::funs(foo)) %>%
        dplyr::mutate(id_str = .data$id_str %>% unlist(recursive = F)) %>% # will be needed from time to time
        dplyr::mutate_(.dots = setNames(list(lazyeval::interp(~ b, b = as.name("_id"))), "underscore_id")) %>% # NSE
        dplyr::select_(lazyeval::interp(~ -a, a = as.name("_id"))) %>% # try harder to understand them!
        dplyr::group_by(.data$id_str) %>% dplyr::slice(1) %>% dplyr::group_by() %>% # there are some "duplicate" tweets w/ slightly different underscore_id
        dplyr::select(-.data$in_reply_to_status_id, -.data$in_reply_to_user_id, -.data$quoted_status_id) # to prevent parsing failures when subsequent tibbles are reread in reread.R
    if ("limit" %in% allnames) {
        df0.tweets %<>%
            dplyr::filter(.data$limit %>% purrr::map_lgl(identical, as.character(NA))) %>% # we want tweets w/ NA limit
            dplyr::select(-.data$limit)
    }
    n.tweets <- df0.tweets %>% nrow
    ## split, mutate those mutatable, & join those unmutatable back
    df1.tweets <- df0.tweets %>%
        ## those mutatable/unlistable
        dplyr::select(-.data$entities, -.data$user, -.data$place, -.data$retweeted_status, -.data$display_text_range, -.data$extended_tweet, -.data$extended_entities, -.data$quoted_status, -.data$geo, -.data$coordinates) %>%
        dplyr::mutate_all(dplyr::funs(unlist, .args = list(recursive = F))) %>%
        dplyr::mutate(time = anytime::anytime(as.double(.data$timestamp_ms) / 1000.0)) %>%
        ## user information
        dplyr::left_join(
                   df0.tweets %>%
                   dplyr::transmute(id_str = .data$id_str, user_list = .data$user %>% purrr::map(quux)) %>% # convert user lists, very important!
                   tidyr::unnest(.data$user_list),
                   by = "id_str"
               ) %>%
        ## retweeted status
        dplyr::left_join(
                   df0.tweets %>% dplyr::transmute(id_str = .data$id_str, .data$retweeted_status),
                   by = "id_str"
               ) %>%
        dplyr::sample_n(n.tweets) %>% # to show rank() is correct
        dplyr::mutate(
                   is_retweet = !is.na(.data$retweeted_status), # consistent w/ retweets
                   original_id_str = .data$retweeted_status %>% purrr::map("id_str") %>% foo %>% unlist, # id of orig. tweet if is retweet, NA if is orig. tweet
                   is_reply =
                       !is.na(.data$in_reply_to_screen_name) &
                       !is.na(.data$in_reply_to_user_id_str) &
                       !is.na(.data$in_reply_to_status_id_str),
                   count = rank(.data$time) %>% as.integer
               )
    df2.tweets <- df1.tweets %>%
        dplyr::left_join(
                   dplyr::filter(df1.tweets, !.data$is_retweet) %>%
                   dplyr::arrange(.data$time) %>%
                   dplyr::transmute(
                              id_str = .data$id_str,
                              interarrival = round(as.numeric(.data$time - lag(.data$time), unit = "secs"), 3L),
                              survival = 1.0 - (rank(.data$interarrival) - 0.5) / sum(!.data$is_retweet),
                              count_orig = rank(.data$time) %>% as.integer
                          ),
                   by = "id_str"
               )
    df3.tweets <- df2.tweets %>%
        dplyr::left_join(
                   dplyr::count(df2.tweets, .data$original_id_str), by = c("id_str" = "original_id_str")
               ) %>%
        ## work w/ counts to prevent parsing failure
        dplyr::mutate(
                   favorite_count = .data$favorite_count %>% as.integer,
                   user_followers_count = .data$user_followers_count %>% as.integer,
                   user_statuses_count = .data$user_statuses_count %>% as.integer,
                   user_friends_count = .data$user_friends_count %>% as.integer,
                   user_favourites_count = .data$user_favourites_count %>% as.integer,
                   user_listed_count = .data$user_listed_count %>% as.integer,
                   retweet_count = ifelse(!.data$is_retweet & is.na(.data$n), 0L, .data$n)
               ) %>%
        dplyr::select(-.data$n, -.data$underscore_id, -.data$user_id)
    ## 02) retweets
    df0.retweets <- df3.tweets %>%
        dplyr::filter(.data$is_retweet) %>%
        dplyr::transmute(retweeted_content = purrr::map2(.data$id_str, .data$retweeted_status, function(x, y) c("retweet_id_str" = x, y))) %>%
        magrittr::extract2("retweeted_content") %>% # now a list of same str. as tweets above!
        purrr::map_df(purrr::map, list) %>% dplyr::mutate_all(dplyr::funs(foo)) %>%
        dplyr::select(-.data$in_reply_to_status_id, -.data$in_reply_to_user_id, -.data$quoted_status_id) # to prevent parsing failures
    n.retweets <- df0.retweets %>% nrow
    df0.retypes <- df0.retweets %>% dplyr::mutate_all(dplyr::funs(purrr::map_chr, .args = list(.f = typeof)))
    check.quoted_status <- magrittr::is_in("quoted_status", df0.retypes %>% names)
    ## split, mutate those mutatable, & join those unmutatable back
    df1.retweets <-
        if (check.quoted_status) {
            df0.retweets %>%
                dplyr::select(-.data$entities, -.data$user, -.data$place, -.data$display_text_range, -.data$extended_tweet, -.data$extended_entities, -.data$geo, -.data$coordinates, -.data$quoted_status) %>%
                dplyr::mutate_all(dplyr::funs(unlist, .args = list(recursive = F)))
        } else {
            df0.retweets %>%
                dplyr::select(-.data$entities, -.data$user, -.data$place, -.data$display_text_range, -.data$extended_tweet, -.data$extended_entities, -.data$geo, -.data$coordinates) %>%
                dplyr::mutate_all(dplyr::funs(unlist, .args = list(recursive = F)))
        } %>%
        dplyr::sample_n(.data$n.retweets) # to show rank() is correct
    df1.retweets %<>%
        ## 01) timestamps for retweets - guaranteed to exist
        dplyr::left_join(df3.tweets %>% dplyr::select(.data$id_str, .data$time), by = c("retweet_id_str" = "id_str")) %>%
        dplyr::rename(retweet_time = .data$time) %>%
        ## 02) timestamps & retweet_counts for original tweets - not guaranteed to exist, as some were created b4 we started
        dplyr::select(-.data$retweet_count) %>% # don't know what this original variable stands for
        dplyr::left_join(df3.tweets %>% dplyr::select(.data$id_str, .data$time, .data$retweet_count), by = c("id_str" = "id_str")) %>%
        ## 03) cum count & survival by original tweet
        dplyr::mutate(
                   count = rank(.data$retweet_time) %>% as.integer,
                   diff = as.numeric(.data$retweet_time - .data$time, unit = "secs") # non-relational
               ) %>%
        dplyr::arrange(.data$retweet_time) %>% # relational below, hence sort 1st
        dplyr::group_by(.data$id_str) %>%
        dplyr::mutate(
                   count_indiv = rank(.data$retweet_time) %>% as.integer, # as.integer to prevent parsing failures
                   favorite_count = .data$favorite_count %>% as.integer, # as.integer to prevent parsing failures
                   interarrival = as.numeric(.data$retweet_time - lag(.data$retweet_time), unit = "secs") %>% round(3),
                   interarrival = ifelse(.data$count_indiv == 1, .data$diff, .data$interarrival),
                   survival_indiv = .data$interarrival %>% qux
               ) %>%
        dplyr::group_by() %>% # same as ungroup()
        ## 04) survival as a whole
        dplyr::mutate(survival = .data$interarrival %>% qux)
    ## 03) original tweets (w/ their retweets)
    t_0 <- df3.tweets$time %>% min %>% lubridate::floor_date(unit = "minutes") # datetime object
    t_inf <- df3.tweets$time %>% max %>% lubridate::ceiling_date(unit = "minutes") %>% magrittr::subtract(t_0) %>% as.numeric(units = "secs") # scalar in secs
    df0.original <- df3.tweets %>%
        dplyr::filter(!.data$is_retweet) %>%
        dplyr::select(.data$time, .data$id_str, .data$user_followers_count, .data$retweet_count) %>%
        dplyr::left_join(
                   df1.retweets %>%
                   dplyr::select(.data$time, .data$id_str, .data$retweet_time),
                   by = c("time", "id_str")
               ) %>%
        dplyr::mutate(t_i = .data$time %>% subtract(t_0) %>% as.numeric(units = "secs"),
                      t_ij = .data$retweet_time %>% subtract(t_0) %>% as.numeric(units = "secs"),
                      t_relative = .data$t_ij - .data$t_i,
                      retweeted = .data$retweet_count > 0L,
                      log_followers_count = log(.data$user_followers_count+1L),
                      log_retweet_count = log(.data$retweet_count+1L)) %>%
        dplyr::select(-.data$time, -.data$retweet_time)
    ## nest t_ij & mutate
    df1.original <- df0.original %>%
        tidyr::nest(.data$t_ij, .data$t_relative) %>%
        dplyr::mutate(t_ij = .data$data %>% purrr::map("t_ij"),
                      t_relative = .data$data %>% purrr::map("t_relative")) %>%
        dplyr::select(-.data$data)
    ## If vec. contains NA, it'll be a single NA
    n.original <- df1.original %>% nrow
    ## checks
    df0.original.a <- df0.original %>% dplyr::filter(is.na(.data$t_ij)) %>% dplyr::arrange(.data$id_str, .data$t_i) # those w/o RTS
    df0.original.b <- df0.original %>% dplyr::filter(!is.na(.data$t_ij)) %>% dplyr::arrange(.data$id_str, .data$t_i, .data$t_ij) # those w/ RTS
    df1.original <-
        dplyr::bind_rows( # restructure
                   df0.original.a %>%
                   dplyr::mutate(
                              t_ij = .data$t_ij %>%
                                  purrr::map(as.numeric),
                              t_relative = .data$t_relative %>%
                                  purrr::map(as.numeric)
                          ),
                   df0.original.b %>%
                   tidyr::nest(.data$t_ij, .data$t_relative) %>%
                   dplyr::mutate(
                              t_ij = .data$data %>%
                                  purrr::map("t_ij"),
                              t_relative = .data$data %>%
                                  purrr::map("t_relative")
                          ) %>%
                   dplyr::select(-.data$data)
               )
    df1.original %>%
        dplyr::mutate(a = .data$t_ij %>% purrr::map_int(length),
                      b = .data$t_ij %>% purrr::map_lgl(identical, as.numeric(NA)),
                      c = ifelse(.data$a == 1L & .data$b, 0L, .data$a)) %>%
        dplyr::count(.data$c == .data$retweet_count) # expect all TRUE - all retweet_counts remain unchanged
    ## save and return
    df0.range <- tibble::data_frame(t_0 = t_0, t_inf = t_inf)
    df0.range %>%
        readr::write_csv(glue("data/{date_str}_range.csv"))
    df0.original.a %>%
        readr::write_csv(glue("data/{date_str}_original_a.csv"))
    df0.original.b %>%
        readr::write_csv(glue("data/{date_str}_original_b.csv"))
    l <- list(
        df0.range = df0.range,
        df0.original.a = df0.original.a,
        df0.original.b = df0.original.b,
        df1.original = df1.original
    )
    invisible(l)
}

#' Re-read files if exist, extract if not
#'
#' @param date_str a character string representing the date.
#' @importFrom readr col_datetime col_double col_character
#' @importFrom dplyr bind_rows mutate select
#' @importFrom purrr map
#' @export
reread <- function(date_str) {
    all_extract_exists <- all(
        file.exists(file_name("range", date_str)),
        file.exists(file_name("original_a", date_str)),
        file.exists(file_name("original_b", date_str))
    )
    if (!all_extract_exists) {
        l <- extraction(date_str)
    } else {
        df0.range <- read_data_csv(
            "range",
            date_str,
            t_0 = readr::col_datetime(),
            t_inf = readr::col_double()
        )
        df0.original.a <- read_data_csv(
            "original_a",
            date_str,
            id_str = readr::col_character(),
            t_ij = readr::col_double(),
            t_relative = readr::col_double(),
            log_retweet_count = readr::col_double()
        )
        df0.original.b <- read_data_csv(
            "original_b",
            date_str,
            id_str = readr::col_character(),
            t_ij = readr::col_double(),
            t_relative = readr::col_double()
        )
        dfa <- dplyr::mutate(df0.original.a, t_ij = purrr::map(.data$t_ij, as.numeric), t_relative = purrr::map(.data$t_relative, as.numeric))
        dfb <- dplyr::select(dplyr::mutate(tidyr::nest(df0.original.b, .data$t_ij, .data$t_relative), t_ij = purrr::map(.data$data, "t_ij"), t_relative = purrr::map(.data$data, "t_relative")), -.data$data)
        df1.original <- dplyr::bind_rows(dfa, dfb)
        l <- list(
            df0.range = df0.range,
            df0.original.a = df0.original.a,
            df0.original.b = df0.original.b,
            df1.original = df1.original
        )
    }
    invisible(l)
}
