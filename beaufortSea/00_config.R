## NOTE: Requires sbtools v1.3.0 or newer
library(sbtools)
source("src/sb_utilities.R")

# Following DSP manual, instructions for using ScienceBase programatically, 
# as of 2/21/2024
## NOTE: Requires sbtools v1.3.0 or newer
initialize_and_cache_sb(
  sb_username = rstudioapi::showPrompt(title = "Username", message = "Username"), # Change to your ScienceBase login
  renviron_file = ".Renviron",
  override_login = FALSE
) 
