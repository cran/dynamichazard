#' @title Hard Drive Failures
#'
#' @description
#' A data set containing hard drive failures data from Backblaze in the
#' start-stop format used in the \code{survival} package.
#'
#' @details
#'
#' Details about the the SMART attributes can be found on
#' \url{https://en.wikipedia.org/wiki/S.M.A.R.T.}. As stated in the original
#' source
#'
#' "Reported stats for the same SMART stat can vary in meaning based on the
#' drive manufacturer and the drive model. Make sure you are comparing
#' apples-to-apples as drive manufacturers don't generally disclose what their
#' specific numbers mean."
#'
#' There are some notes on \url{https://en.wikipedia.org/wiki/S.M.A.R.T.}
#' regarding which attributes that have vendor specific raw value. Further,
#'
#' "The values in the files are the values reported by the drives. Sometimes,
#' those values are out of whack. For example, in a few cases the RAW value of
#' SMART 9 (Drive life in hours) reported a value that would make a drive 10+
#' years old, which was not possible. In other words, it's a good idea to have
#' bounds checks when you process the data."
#'
#' See this github page for the processing steps
#' \url{https://github.com/boennecd/backblaze_survival_analysis_prep}.
#'
#'
#' @format A \code{data.frame} with the following columns:
#' \describe{
#'   \item{serial_number}{Serial number for the hard disk which the row belongs
#'   to.}
#'   \item{model}{hard disk model.}
#'   \item{manufacturer}{manufacturer of the hard disk model.}
#'   \item{tstart,tstop}{start and stop times on the SMART 9 attribute scale.}
#'   \item{fails}{1 if the hard disk fails at \code{tstop}.}
#'   \item{size_tb}{hard disk size in terabytes.}
#'   \item{smart_x}{the raw SMART attribute x value. E.g., \code{smart_12} is the power cycle count.}
#'   \item{smart_x_bin}{1 if the SMART attribute x value is non-zero.}
#'   \item{..._cumsum}{cumulative sum of the prefix \code{...}.}
#'   \item{n_fails}{number of failures in the original data. Hard disk should
#'   only fail once but this is not the case in the raw data.}
#'   \item{n_records}{number of records in the original source.}
#'   \item{min_date,max_date}{first and last date in the original source.}
#'   \item{min_hours,max_hours}{smallest and largest value of the SMART 9
#'   attribute in the original source.}
#' }
#'
#' @source
#' Raw data from \url{https://www.backblaze.com/b2/hard-drive-test-data.html}.
#' Data have been processed to get a start-stop \code{data.frame} format.
"hds"
