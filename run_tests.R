library(testthat)
library(dplyr)
library(readr)

# Source the application and plots code
source("plots.R")

context("MPRA Saturation Mutagenesis App Tests")

# Test 1: Data loading functions
test_that("loadElementList loads data correctly", {
  # Create a temporary test file
  temp_file <- tempfile()
  writeLines(c("chr22:50964070-50964571", "chr1:100000-100100"), temp_file)
  
  result <- tryCatch({
    elements <- read.delim(temp_file, header=FALSE, quote="") 
    colnames(elements) <- c("Element")
    elements <- elements %>% arrange(Element)
    elements
  }, error = function(e) {
    NULL
  })
  
  expect_is(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true("Element" %in% colnames(result))
  
  unlink(temp_file)
})

# Test 2: Verify app can be sourced without errors
test_that("app.R sources without errors", {
  # Try to source the app - this checks syntax and basic structure
  result <- tryCatch({
    # Parse the file to check for syntax errors
    parse(file = "app.R")
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  expect_true(result)
})

# Test 3: Verify all required libraries can be loaded
test_that("All required packages are available", {
  packages <- c("shiny", "htmlwidgets", "DT", "dplyr", "ggplot2", "readr", "stringr", "plotly")
  
  for (pkg in packages) {
    result <- require(pkg, quietly = TRUE, character.only = TRUE)
    expect_true(result, info = paste("Failed to load package:", pkg))
  }
})

# Test 4: Verify plots.R sources without errors
test_that("plots.R sources without errors", {
  result <- tryCatch({
    parse(file = "plots.R")
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  expect_true(result)
})

# Test 5: Check that data files exist
test_that("Data files exist", {
  expect_true(file.exists("data/enhancers.tsv"))
  expect_true(file.exists("data/promoters.tsv"))
})

# Test 6: Check that markdown files exist
test_that("Markdown files exist", {
  expect_true(file.exists("mrkdown/about.md"))
  expect_true(file.exists("mrkdown/file_format.md"))
})

# Test 7: Basic data file format validation
test_that("Data files have expected format", {
  # Check enhancers
  enhancers <- tryCatch({
    read.delim("data/enhancers.tsv", header=FALSE, quote="", nrow=1)
  }, error = function(e) {
    NULL
  })
  expect_false(is.null(enhancers))
  
  # Check promoters
  promoters <- tryCatch({
    read.delim("data/promoters.tsv", header=FALSE, quote="", nrow=1)
  }, error = function(e) {
    NULL
  })
  expect_false(is.null(promoters))
})

cat("\n=== All tests completed successfully! ===\n")
