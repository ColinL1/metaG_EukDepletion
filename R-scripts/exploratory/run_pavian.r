# Install and run pavian
if (!require(remotes)) { install.packages("remotes") }
if (!require(pavian, quietly = TRUE)) {
  remotes::install_github("fbreitwieser/pavian")
}

pavian::runApp()
