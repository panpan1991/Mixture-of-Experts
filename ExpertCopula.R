h3=function(x){
  1/4*(1/0.05*dnorm((x-0.15)/0.05)+1/0.2*dnorm((x-0.6)/0.2))
}
#generate data
#x=runif(300)
x=seq(0,  1, length.out = 100)
n=length(x)
y=rnorm(n,h3(x), 0.5)
hist(y)
#Scatter plot with the real line
plot(x,y)
par(new=TRUE)
plot(x,h3(x), type = 'l')

data=data.frame(x,y)

write.table(data, 'h3data', row.names = TRUE, col.names = TRUE)









