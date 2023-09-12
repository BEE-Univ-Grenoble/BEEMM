#' @import readr
#' @import stringr
#' @import dplyr
#' @import tibble
#' @author Eric Coissac
NULL

#' Read DNA sequences from a file
#'
#' @author Eric Coissac
#' @md
#'
#' read_fasta() reads up to n_max DNA sequences from a fasta file.
#'
#' @param file Either a path to a file, a connection, or literal data (either a
#'   single string or a raw vector). Files ending in .gz, .bz2, .xz, or .zip
#'   will be automatically uncompressed. Files starting with  ⁠http://⁠,
#'   ⁠https://⁠, ⁠ftp://⁠, or  ⁠ftps://⁠ will be automatically downloaded.
#'   Remote gz files can also be automatically downloaded and decompressed.
#'   Literal data is most useful for examples and tests. To be recognised as
#'   literal data, the input must be either wrapped with I(), be a string
#'   containing at least one new line, or be a vector containing at least one
#'   string with a new line. Using a value of clipboard() will read from the
#'   system clipboard.
#' @param progress Display a progress bar? By default it will only display in an
#'   interactive session and not while knitting a document. The automatic
#'   progress bar can be disabled by setting option readr.show_progress to
#'   FALSE.
#' @param num_threads The number of processing threads to use for initial
#'   parsing.
#'
#' @return a tibble with two columns
#'   - `id` : the sequence identifier (character)
#'   - `sequence` : the DNA sequence (character)
#'
#' @export
#'
#' @examples
read_fasta <- function(file, progress = FALSE,
                       num_threads = readr_threads()) {
  read_lines(file,
    skip_empty_rows = TRUE,
    lazy = FALSE,
    progress = progress,
    num_threads = num_threads
  ) %>%
    tibble(line = .) %>%
    mutate(
      beginning = str_detect(line, pattern = "^>"),
      number = cumsum(beginning)
    ) %>%
    group_by(number, beginning) %>%
    summarise(lines = str_c(line, collapse = "")) %>%
    ungroup() %>%
    mutate(attribut = ifelse(beginning, "id", "sequence")) %>%
    select(-beginning) %>%
    pivot_wider(
      id_cols = number,
      names_from = attribut,
      values_from = lines
    ) %>%
    select(id, sequence) %>%
    mutate(id = str_extract(id, pattern = "(?<=^>)[:graph:]+"))
}
