
ali <- 330
david <- 8128
ben <- 1037

set.seed(ali+ben+david)

print(c(
	"AG", "FA", "DE", 
	sample(c("BB", "JD"))
))
