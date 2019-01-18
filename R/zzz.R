.onAttach <- function(libname, pkgname) {
  setHook(packageEvent("data.table", "attach"), function(...) {
    packageStartupMessage(
      "You have loaded data.table after exSTRa. As exSTRa replaces copy() with a generic,\n",
      "this may cause issues when copy()ing exSTRa objects.")
  })
}