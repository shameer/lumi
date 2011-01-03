setGeneric("methylated", function(object) standardGeneric("methylated"))
setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
setGeneric("getHistory", function(object) standardGeneric("getHistory"))
setGeneric("MAplot", function(object, ...) standardGeneric("MAplot"))

# if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
# if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
# if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
# if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
# if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
# if (is.null(getGeneric("MAplot"))) setGeneric("MAplot", function(object, ...) standardGeneric("MAplot"))
# if (is.null(getGeneric("plotCDF"))) setGeneric("plotCDF", function(x, ...) standardGeneric("plotCDF"))

if (is.null(getGeneric("beadNum"))) setGeneric("beadNum", function(object) standardGeneric("beadNum"))
if (is.null(getGeneric("beadNum<-"))) setGeneric("beadNum<-", function(object, value) standardGeneric("beadNum<-"))

if (is.null(getGeneric("detection"))) setGeneric("detection", function(object) standardGeneric("detection"))
if (is.null(getGeneric("detection<-"))) setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))

if (is.null(getGeneric("controlData"))) setGeneric("controlData", function(object) standardGeneric("controlData"))
if (is.null(getGeneric("controlData<-"))) setGeneric("controlData<-", function(object, value) standardGeneric("controlData<-"))

if (is.null(getGeneric("summary"))) setGeneric("summary", function(object, ...) standardGeneric("summary"))
if (is.null(getGeneric("show"))) setGeneric("show", function(object) standardGeneric("show"))
if (is.null(getGeneric("combine"))) setGeneric("combine", function(x, y, ...) standardGeneric("combine"))
