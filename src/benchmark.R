result <- microbenchmark(
    Refactored = source("src/AlienDetective.R"),
    MGielen = source("src/AlienDetective_MGielen.R"),
    times = 1
  )

print(result)
