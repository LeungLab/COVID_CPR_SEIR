

data <- read.csv("MainModel.csv")
data0 <- data[data$true == 0, ]
data1 <- data[data$true == 1, ]

# define colors used in plots
DarkGreen = "#3A5431"
CustomGrey = "#858585"

ppl.wanting.tests <- 5000


####### DEFINE FUNCTION ##########

daily.fun <- function(prop.true.pos,lgnd) {
  n0 <-
    round(ppl.wanting.tests * (1 - prop.true.pos))# number of 0s that want to get tested
  n1 <-
    round(ppl.wanting.tests * prop.true.pos)  # number of 1s that want to get tested
  
  
  proportions <- function(Ntests) {
    # set up patient pool of people who all want testing
    patients0 <- data0[sample(nrow(data0), n0, replace = TRUE), ]
    patients1 <- data1[sample(nrow(data1), n1, replace = TRUE), ]
    samp <- rbind(patients0, patients1)
    
    samp <- samp[order(-samp$pred_RF), ]
    
    # how many 1s do we test if we sample using the model?
    sum(samp[1:Ntests, ]$true == 1)
    prop.preferred.testing <- sum(samp[1:Ntests, ]$true == 1) / n1
    
    # how many 1s do we test if we sample randomly
    sum(samp[sample(nrow(samp), Ntests), ]$true == 1)
    prop.random.testing <-
      sum(samp[sample(nrow(samp), Ntests), ]$true == 1) / n1
    
    return(c(prop.preferred.testing, prop.random.testing))
  }
  
  reps <- function(Ntests) {
    rowMeans(replicate(1000, proportions(Ntests))) # should be 1000
  }
  
  Ntests.range <- round(seq(0, 1, .02) * (n0 + n1))
  results <- sapply(Ntests.range, reps)
  
  ### plot results ###
  {
    par(mar = c(3,3,3,3)) # margin: bottom, left, top, right
    cex.val = 1.25
    plot(
      100 * Ntests.range / (n0 + n1),
      results[1, ] * 100,
      type = "l",
      col = DarkGreen,
      lwd = 2,
      xlab="", ylab="",
      ylim = c(0, 100),
      cex.main = cex.val,
      cex.lab = cex.val,
      cex.axis = cex.val
    )
    lines(100 * Ntests.range / (n0 + n1),
          results[2, ] * 100,
          col = CustomGrey,
          lwd = 2)
    
    par(new = T)
    plot(
      100 * Ntests.range / (n0 + n1),
      results[1, ] / (results[2, ]),
      type = "l",
      lty = "dotted",
      lwd = 2,
      axes = F,
      xlab = NA,
      ylab = NA,
      ylim = c(0, 5),
      cex.lab = 1.25,
      cex.axis = 1.25
    )
    axis(side = 4, cex.axis = 1.25)
  }
}

###### make plots #####
filename = "Individual_day_results.pdf"
pdf(file=filename,width=7.5, height=11*.6)

# specify layout of multiple subplots
par(oma=c(2,2,1,2)) # bottom, left, top, right
m <- matrix(c(1,1,1,1,
              2,2,3,3,
              4,4,5,5),
            3,4, byrow = TRUE)
layout(m, heights = c(0.07,0.3,0.3))

### put legend at top ###
{
  par(mar = c(.2,1,1,1))# bottom, left, top, right
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x="bottomright", inset=c(0,0),
        legend=c("prioritized testing", "indiscriminate testing", "fold change"),
        text.width=c(0.154,0.15,0.1),
        lwd = c(2, 2, 2),
        bty = "n",
        col = c(DarkGreen, CustomGrey, "black"),
        lty = c("solid", "solid", "dotted"),
        cex = 1.5, horiz=TRUE)
}

{
### 5% ###
daily.fun(0.05,lgnd=TRUE)
mtext("A", side=3, adj=0, line=.5, cex=1.2); 
param.txt <- bquote(paste('q = ', .(0.05)))
mtext(param.txt, side=3, adj=1, line=.5, cex=1.2)

### 25% ###
daily.fun(0.25,lgnd=FALSE)
mtext("B", side=3, adj=0, line=.5, cex=1.2); 
param.txt <- bquote(paste('q = ', .(0.25)))
mtext(param.txt, side=3, adj=1, line=.5, cex=1.2)

### 50% ###
daily.fun(0.5,lgnd=FALSE)
mtext("C", side=3, adj=0, line=.5, cex=1.2); 
param.txt <- bquote(paste('q = ', .(0.50)))
mtext(param.txt, side=3, adj=1, line=.5, cex=1.2)

### 75% ###
daily.fun(0.75,lgnd=FALSE)
mtext("D", side=3, adj=0, line=.5, cex=1.2); 
param.txt <- bquote(paste('q = ', .(0.75)))
mtext(param.txt, side=3, adj=1, line=.5, cex=1.2)
}

mtext("% of test eligible people tested",side=1,line=0,outer=TRUE,cex=1.2)
mtext("% of SARS-CoV-2 positive people tested",side=2,line=0,adj=0.35,outer=TRUE,cex=1.2,las=0)
mtext("fold change",side = 4,line = 0,adj=0.44,outer=TRUE,cex=1.2)
dev.off()
