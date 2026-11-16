# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# The thread pool's threads run code in ecodive's DLL, so they are stopped
# before the DLL is unloaded (as by unloadNamespace() or devtools::load_all()).
.onUnload <- function (libpath) {
  .Call(C_stop_thread_pool)
  library.dynam.unload('ecodive', libpath)
}
