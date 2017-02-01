#Click the Source button to not echo the commands.
finish <- 1
while(finish == 1){
  #Getting Started
  if(7+2==9){
    print("Testing")
  }else{
    break}
  x <- 7+2
  if(x^2==81){
    print("Systems?")
    print("   Check!")
  }else{
    break}
  
  #Basic Commands
  y <- sqrt(1.21)
  
  #is pi exact?
  if(pi == 3.141593){
    #It isn't
    print("Captain Wilco")
  }else{
    #It is more exact than this
    print("Captain Solo")
  }
  print("Engines?")
  
  #Exercise 1
  s <- sin(pi/3) + log(3)
  if(s^4/2 < exp(s)){
    print("   Roger")
    print("")
    print("Exercise 1: Conclusion")
    print("s^4/2 < exp(s)")
    print("")
  }else{
    print("   Yes, See Three!")
    print("")
    print("Exercise 1: Conclusion")
    print("s^4/2 >= exp(s)")
    print("")
  }
  
  a <- 1+2i
  b <- 1i
  #Does R tolerate errors in the middle of scripts?
  #No!
  #c <- i
  
  if(a/b == 2-1i){
    print("T Minus")
  }#So it has some flexibility with placement of }?
  else{
    break
  }
  
  if(Conj(a) == 1-2i){
    print("10")
  }else{
    break
  }
  
  #Vectors
  #Space tolerant?
  x<-c(      9,8 ,7, 6)
  
  #Entry recycling is kind of bizarre, but good to know.
  cat(x / (1:8) + c(10,200), "\n")
  print("   ...Houston forgot how to count?")
  
  #Checking how many entries are checked in conditional
  #Warnings are tolerated, but printed at end of script
  print("S")
  if(x*c(TRUE, FALSE)){
    print("tut")
  }
  print("ter?")
  
  #Display a long vector
  print(2^seq(10,0, by = -.2))
  print("   ...Houston...?")
  
  #Exercise 2
  t <- seq(1,10)
  print("")
  print("Exercise 2: Conclusion")
  cat('     t cos(t) = ', t*cos(t), '\n     t^2 exp(-t) = ', t^2 * exp(-t), '\n')
  print("")
  
  
  #Exercise 3
  print("Exercise 3: Conclusion")
  for (i in seq(0,20)){
    n <- seq(0,i)
    r <- .5
    x <- r^n
    cat(i,' ', sum(x),' ', sum(x)-(1), '\n')
  }
  cat('inf ',1/(1-r), ' 0\n')
  print("")
  
  #Can we just print matrices?
  A <- matrix(seq(1,9), nrow = 3, ncol = 3)
  print(A <- cbind(rbind(t(A),c(1,1,1)),c(1,2,3,4),deparse.level=0))
  #Note that traversal is top to bottom first
  #We also can do assignment in the middle of a command!
  print(eigen(A))
  
  #Exercise 4
  print("Exercise 4")
  print(A <- t(matrix(c(1,-2,3,5, -2, 1, 1, 1, -1, 3, 1, -4, 1, 1, 1, 7), nrow = 4, ncol = 4)))
  RightVal <- c(-4,5, 13, -4)
  print(solve(A,RightVal))
  
  #Evidently forgot to push on laptop
  #Basic Plotting
  print("We plot in the (x,y)-plane via plot(x,y)")
  x <- seq(0, 2*pi, by=.1)
  y <- sin(x)
  plot(x,y, type="o")
  #Note that changing plot parameters
  #only affects future plots!
  #Parameters do stay in place between
  #runs.
  
  #Exercise 5
  par(pty="s", bg = "black", col.lab = "white", col.main = "white", col.axis = "white")
  theta=(0:100)*2*pi/100;
  x <- cos(theta)
  y <- sin(theta) 
  coloring <- rgb(x*.5+.5,y*.5+.5,0*theta)
  plot(x, y, main="Exercise 5", xlab="X", ylab="y", col = coloring)
  grid(col = "white")
  
  #Cleanup -> Reset parameters
  par(pty="m")
  
  finish <- 0
}
if(finish == 1){
  print("Captain, there's been an error. Launch aborted")
}else{#Note, if the } is not adjacent to the else, 
  #R does not detect as the same if-else
  print("Liftoff!")
  print("   We have Liftoff!")}


