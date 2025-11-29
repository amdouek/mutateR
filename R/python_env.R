#' Install the mutateR Python environment
#'
#' Sets up a specific Conda environment for mutateR containing
#' TensorFlow and necessary dependencies.
#'
#' @param envname Character. Name of the conda environment (default "r-mutater").
#' @param python_version Character. Python version to install (default "3.9").
#' @param fresh Logical. If TRUE, removes existing environment before installing.
#'
#' @export
install_mutater_env <- function(envname = "r-mutater",
                                python_version = "3.9",
                                fresh = FALSE) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required.")
  }

  # 1. Check/Install Miniconda
  tryCatch({
    conda_path <- reticulate::miniconda_path()
    if (!file.exists(conda_path)) {
      message("Miniconda not found. Installing Miniconda via reticulate...")
      reticulate::install_miniconda()
    }
  }, error = function(e) {
    stop("Error checking Miniconda: ", e$message)
  })

  # 2. Check existing environment
  envs <- reticulate::conda_list()
  env_path <- envs$python[envs$name == envname]

  if (length(env_path) > 0) {
    if (fresh) {
      message("Removing existing environment '", envname, "'...")
      reticulate::conda_remove(envname)
    } else {
      message("Environment '", envname, "' already exists. Set fresh=TRUE to reinstall.")
      return(invisible(NULL))
    }
  }

  # 3. Create Environment
  message("Creating conda environment '", envname, "' with Python ", python_version, "...")
  reticulate::conda_create(envname, python_version = python_version)

  # 4. Install Dependencies
  # Note: Pinning numpy<2 is currently recommended for TensorFlow compatibility
  # Using tensorflow-cpu for Windows compatibility.
  # If on Mac M1/M2, users might need to manually install tensorflow-macos.
  pkgs <- c("tensorflow-cpu", "numpy<2", "h5py", "pandas", "scipy")

  message("Installing packages: ", paste(pkgs, collapse = ", "))
  reticulate::conda_install(envname, packages = pkgs, pip = TRUE)

  message("\nInstallation complete! Please restart your R session before activating.")
}

#' Activate the mutateR Python environment
#'
#' Points the current R session to the mutateR Conda environment.
#'
#' @param envname Character. Name of the environment (default "r-mutater").
#'
#' @return Logical TRUE if successful, FALSE otherwise.
#' @export
activate_mutater_env <- function(envname = "r-mutater") {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)

  # 1. Check if Reticulate is ALREADY initialized
  if (reticulate::py_available()) {
    config <- reticulate::py_config()

    # Robust check: does the path contain our environment name?
    # (reticulate config$envname is often the full path, not just the name)
    is_correct <- config$envname == envname || grepl(envname, config$python, fixed = TRUE)

    if (!is_correct) {
      warning("Reticulate is already initialized to a different environment: '", config$envname, "'.\n",
              "You must RESTART R to switch to '", envname, "'.")
      return(FALSE)
    }
    return(TRUE)
  }

  # 2. Check RETICULATE_PYTHON collision
  sys_ret_py <- Sys.getenv("RETICULATE_PYTHON")

  if (sys_ret_py != "") {
    if (!grepl(envname, sys_ret_py, fixed = TRUE)) {
      message(sprintf(
        "NOTE: RETICULATE_PYTHON is set to '%s'. This may conflict with mutateR.\n",
        sys_ret_py
      ))
      # We cannot reliably unset it if reticulate already loaded,
      # so we warn the user instead of silently failing.
    }
  }

  # 3. Attempt Activation
  tryCatch({
    reticulate::use_condaenv(envname, required = TRUE)

    # 4. Verify
    config <- reticulate::py_config()

    if (!grepl(envname, config$python, fixed = TRUE)) {
      # Fallback check
      if(config$envname != envname) {
        warning("Failed to activate '", envname, "'. Active python: ", config$python)
        return(FALSE)
      }
    }

    # 5. Check TensorFlow
    if (!reticulate::py_module_available("tensorflow")) {
      warning("Environment activated, but TensorFlow not found.")
      return(FALSE)
    }

    return(TRUE)

  }, error = function(e) {
    warning("Activation error: ", e$message)
    return(FALSE)
  })
}

#' Check mutateR environment status
#'
#' Diagnostics for the Python setup.
#'
#' @export
check_mutater_env <- function(envname = "r-mutater") {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate missing.")

  # Check 1: Conda List
  envs <- reticulate::conda_list()
  exists <- envname %in% envs$name

  if (!exists) {
    # If not found by name, check if it's currently active (sometimes conda_list misses paths)
    if (reticulate::py_available()) {
      conf <- reticulate::py_config()
      if (grepl(envname, conf$python, fixed=TRUE)) exists <- TRUE
    }
  }

  if (!exists) {
    return(FALSE)
  }

  # Check 2: Active Session
  # We return TRUE if the env exists on disk AND (if python is initialized) it is the correct one.
  if (reticulate::py_available()) {
    conf <- reticulate::py_config()
    # Check if the path contains the environment name
    is_correct <- grepl(envname, conf$python, fixed = TRUE) || (conf$envname == envname)
    return(is_correct)
  }

  # If Python isn't running yet, but the env exists, we return TRUE
  # (assuming activation will happen later)
  return(TRUE)
}
