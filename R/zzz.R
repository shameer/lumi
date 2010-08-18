.First.lib <- function(libname, pkgname) {
    
    if(.Platform$OS.type == "windows" && interactive()
        && .Platform$GUI ==  "Rgui") {
        addVigs2WinMenu("lumi")
    }

}
