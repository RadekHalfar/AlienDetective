result <- microbenchmark::microbenchmark(
    #Refactored = source("src/AlienDetective.R"),
    Refactored_v2 = source("src/AlienDetective_v2.R"),
    MGielen = source("src/AlienDetective_MGielen.R"),
    times = 3
  )
  
print(result)
