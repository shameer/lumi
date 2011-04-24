setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setGeneric("pairs", function(x, ...) standardGeneric("pairs"))
setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))
setGeneric("density", function(x, ...) standardGeneric("density"))
setGeneric("hist", function(x, ...) standardGeneric("hist"))
setGeneric("MAplot", function(object, ...) standardGeneric("MAplot"))

setGeneric("beadNum", function(object) standardGeneric("beadNum"))
setGeneric("beadNum<-", function(object, value) standardGeneric("beadNum<-"))

setGeneric("controlData", function(object) standardGeneric("controlData"))
setGeneric("controlData<-", function(object, value) standardGeneric("controlData<-"))

setGeneric("detection", function(object) standardGeneric("detection"))
setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))

setGeneric("getHistory", function(object) standardGeneric("getHistory"))


if (is.null(getGeneric("summary"))) setGeneric("summary", function(object, ...) standardGeneric("summary"))
if (is.null(getGeneric("show"))) setGeneric("show", function(object) standardGeneric("show"))
if (is.null(getGeneric("combine"))) setGeneric("combine", function(x, y, ...) standardGeneric("combine"))
