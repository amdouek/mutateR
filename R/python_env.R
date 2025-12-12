#' Install the mutateR Python environment
#'
#' Sets up a specific Conda environment for mutateR
#' containing TensorFlow, primer3-py, and necessary dependencies.
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

  # 3. Create environment
  message("Creating conda environment '", envname, "' with Python ", python_version, "...")
  reticulate::conda_create(envname, python_version = python_version)

  # 4. Install dependencies
  # Note: Pinning numpy<2 is currently recommended for TensorFlow compatibility
  pkgs <- c("tensorflow-cpu", "numpy<2", "h5py", "pandas", "scipy", "primer3-py", "rs3") # We currently don't really need gpu for what we do here, so tensorflow-cpu will suffice

  message("Installing packages: ", paste(pkgs, collapse = ", "))
  reticulate::conda_install(envname, packages = pkgs, pip = TRUE)

  message("\nInstallation complete! Please restart your R session before activating.")
}

#' Activate the mutateR Python environment
#'
#' Points the current R session to the mutateR Conda environment.
#' Automatically handles RETICULATE_PYTHON conflicts.
#'
#' @param envname Character. Name of the environment (default "r-mutater").
#'
#' @return Logical TRUE if successful, FALSE otherwise.
#' @export
activate_mutater_env <- function(envname = "r-mutater") {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)

  # 1. Check if reticulate is already initialised
  if (reticulate::py_available()) {
    config <- reticulate::py_config()

    # Safe string extraction
    curr_env <- if (is.null(config$envname) || is.na(config$envname)) "" else config$envname
    curr_py  <- if (is.null(config$python) || is.na(config$python)) "" else config$python

    # Check if we are already in the correct environment
    is_correct <- (curr_env == envname) || grepl(envname, curr_py, fixed = TRUE)

    if (!is_correct) {
      warning("Reticulate is already initialized to a different environment: '", curr_env, "'.\n",
              "Path: ", curr_py, "\n",
              "You must RESTART R to switch to '", envname, "'.")
      return(FALSE)
    }
    return(TRUE)
  }

  # 2. Check and Fix RETICULATE_PYTHON collision
  sys_ret_py <- Sys.getenv("RETICULATE_PYTHON")
  if (sys_ret_py != "") {
    if (!grepl(envname, sys_ret_py, fixed = TRUE)) {
      message("Conflict detected: RETICULATE_PYTHON is set to '", sys_ret_py, "'.")
      message("Unsetting RETICULATE_PYTHON for this session to allow mutateR environment activation...")
      Sys.unsetenv("RETICULATE_PYTHON")
    }
  }

  # 3. Attempt env activation
  tryCatch({
    reticulate::use_condaenv(envname, required = TRUE)

    # 4. Verify
    config <- reticulate::py_config()
    curr_env <- if (is.null(config$envname) || is.na(config$envname)) "" else config$envname
    curr_py  <- if (is.null(config$python) || is.na(config$python)) "" else config$python

    if (!grepl(envname, curr_py, fixed = TRUE)) {
      # Fallback check
      if(curr_env != envname) {
        warning("Failed to activate '", envname, "'. Active python: ", curr_py)
        return(FALSE)
      }
    }

    # 5. Check TensorFlow (Primary Dependency)
    if (!reticulate::py_module_available("tensorflow")) {
      warning("Environment activated, but TensorFlow not found. Try running install_mutater_env(fresh=TRUE).")
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

  # Check 1: Conda list
  envs <- reticulate::conda_list()
  exists <- envname %in% envs$name

  if (!exists) {
    if (reticulate::py_available()) {
      conf <- reticulate::py_config()
      curr_py <- if (is.null(conf$python) || is.na(conf$python)) "" else conf$python
      if (grepl(envname, curr_py, fixed=TRUE)) exists <- TRUE
    }
  }

  if (!exists) {
    return(FALSE)
  }

  # Check 2: Active session
  if (reticulate::py_available()) {
    conf <- reticulate::py_config()
    curr_env <- if (is.null(conf$envname) || is.na(conf$envname)) "" else conf$envname
    curr_py  <- if (is.null(conf$python) || is.na(conf$python)) "" else conf$python

    is_correct <- grepl(envname, curr_py, fixed = TRUE) || (curr_env == envname)
    return(is_correct)
  }

  # If Python isn't running yet, but the env exists, return TRUE
  return(TRUE)
}
