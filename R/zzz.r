# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# The thread pool's threads are also stopped when R exits, so that tools
# such as valgrind don't report their memory as possibly lost. The finalizer
# does nothing once .onUnload() has stopped the pool and unloaded the DLL.
.onLoad <- function (libname, pkgname) {
  guard <- new.env(parent = emptyenv())
  guard$unloaded <- FALSE
  assign('.pool_guard', guard, envir = topenv())
  reg.finalizer(guard, onexit = TRUE, function (e) {
    if (!isTRUE(e$unloaded)) .Call(C_stop_thread_pool)
  })
}


# The thread pool's threads run code in ecodive's DLL, so they are stopped
# before the DLL is unloaded (as by unloadNamespace() or devtools::load_all()).
.onUnload <- function (libpath) {
  .Call(C_stop_thread_pool)
  guard <- get0('.pool_guard', envir = topenv(), inherits = FALSE)
  if (is.environment(guard)) guard$unloaded <- TRUE
  library.dynam.unload('ecodive', libpath)
}
