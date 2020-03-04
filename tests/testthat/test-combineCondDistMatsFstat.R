## Tests for .combineCondDistMatsFstat
## library(Rtpca); library(testthat)
context(".combineCondDistMatsFstat")

library(dplyr)
library(Biobase)

c1_df <- tibble(
    rowname = "3",
    key = "2",
    value = 2,
    pair = "2:3"
)

c2_df <- tibble(
    rowname = "3",
    key = "2",
    value = 0.9,
    pair = "2:3"
)

ctrl_df <- tibble(
    pair = "2:3",
    valueC1 = 2,
    valueC2 = 0.9,
    rssC1 = sqrt(2),
    rssC2 = sqrt(0.9),
    rssC1_rssC2 = sqrt(2) - sqrt(0.9),
    min_rssC1_rssC2 = sqrt(0.9),
    f_stat = (sqrt(2) - sqrt(0.9))/sqrt(0.9)
)

expect_equal(
    .combineCondDistMatsFstat(c1_df, c2_df),
    ctrl_df
)
