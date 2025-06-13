result <- microbenchmark::microbenchmark(
    Refactored = source("src/AlienDetective.R"),
    MGielen = source("src/AlienDetective_MGielen.R"),
    times = 2
  )

print(result)
