###################
# IsDate function #
###################

IsDate <- function(mydate, date.format = "%d-%m-%y") {
  tryCatch(!is.na(as.Date(mydate, date.format)),  
           error = function(err) {FALSE})  
}