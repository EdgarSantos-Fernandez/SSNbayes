test_that("Boise example files can be located and the observations read", {
  boise_path <- system.file("extdata", "boise.ssn", package = "SSNbayes")
  obs_path <- system.file("extdata", "boise_obs.RDS", package = "SSNbayes")

  expect_true(nzchar(boise_path))
  expect_true(dir.exists(boise_path))
  expect_true(file.exists(file.path(boise_path, "edges.shp")))
  expect_true(file.exists(file.path(boise_path, "sites.shp")))
  expect_true(file.exists(file.path(boise_path, "preds.shp")))
  expect_true(file.exists(file.path(boise_path, "binaryID.db")))

  expect_true(nzchar(obs_path))
  boise_obs <- readRDS(obs_path)
  expect_s3_class(boise_obs, "data.frame")
  expect_identical(dim(boise_obs), c(3654L, 49L))
  expect_true(all(c(
    "pid", "netID", "date", "locID", "temp", "air_temp",
    "slope", "elev", "drainage", "sin", "cos", "dataset"
  ) %in% names(boise_obs)))
  expect_identical(length(unique(boise_obs$locID)), 42L)
  expect_identical(length(unique(boise_obs$date)), 87L)
})
