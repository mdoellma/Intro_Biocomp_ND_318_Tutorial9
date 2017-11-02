#set up
	#load packages
	#read data file
	
#functions
	#def function_name(argument, argument, argument):
		#do stuff here
		#result_variable=do final stuff here
		#return result_variable
	#def nllike(arguments):
		#"unpack" arguments, assign variables
		#calculate expected value (model equation)
		#result = calculated negative log likelihood
		#return result (the negative log likelihood)
	def nllike(p,obs)
		B0=p[0]
		B1=p[1]
		sigma=p[2]
		expected=B0+B1*obs.x
			#from the data set obs, use the column titled x
			#for t-tests, obs.x will only be 0 or 1 (control and treatment)
		nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
			#convert to the log scale
			
#likelihood ratio tests
	#D = 2*(nll alternative model-nll null model)
	#D ~ X^2, df = difference in parameters between models
	#pval = probability of a "more extreme" D if two models are equal
	#we estimate 3 parameters in t-test and 2 parameters in the null model
	
	#code
		#from scipy.stats import chi2
		#pval=1-chi2.df(x=2*...
		#the rest is posted!
			