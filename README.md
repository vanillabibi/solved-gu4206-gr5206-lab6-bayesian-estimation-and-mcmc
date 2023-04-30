Download Link: https://assignmentchef.com/product/solved-gu4206-gr5206-lab6-bayesian-estimation-and-mcmc
<br>
Make sure that you upload the PDF (or HTML) output after you have knitted the file. The files you upload to the Canvas page should be updated with commands you provide to answer each of the questions below. You can edit this file directly to produce your final solutions.

<h2><strong>Goals</strong></h2>

This lab has two goals. The first goal is to use the <strong><sub>Accept-Reject </sub></strong>algorithm to simulate from a mixture of two normals. The second goal is to utilize Bayesian methods and the the famous <strong><sub>Markov Chain Monte </sub>Carlo </strong>algorithm to estimate the mixture parameter <em><sub>δ</sub></em>.

<h1>Background: (Mixture)</h1>

A mixture distribution is the probability distribution of a random variable that is derived from a collection of other random variables (Wiki). In our case we consider a mixture of two normal distributions. Here we assume that our random variable is governed by the probability density <em><sub>f</sub></em>(<em><sub>x</sub></em>), defined by

<em>f</em>(<em>x</em>) = <em>f</em>(<em>x</em>;<em>µ</em><sub>1</sub><em>,σ</em><sub>1</sub><em>,µ</em><sub>2</sub><em>,σ</em><sub>2</sub><em>,δ</em>)

<em>,</em>

where <em>−∞ &lt; x &lt; ∞ </em>and the parameter space is defined by <em>−∞ &lt; µ</em><sub>1</sub><em>,µ</em><sub>2 </sub><em>&lt; ∞</em>, <em>σ</em><sub>1</sub><em>,σ</em><sub>2 </sub><em>&gt; </em>0, and 0 <em>≤ δ ≤ </em>1.

The <strong>mixture parameter </strong><em><sub>δ </sub></em>governs how much mass gets placed on the first distribution <em><sub>f</sub></em>(<em><sub>x</sub></em>;<em><sub>µ</sub></em><sub>1</sub><em><sub>,σ</sub></em><sub>1</sub>) and the complement of <em><sub>δ </sub></em>governs how much mass gets placed on the other distribution <em><sub>f</sub></em><sub>2</sub>(<em><sub>x</sub></em>;<em><sub>µ</sub></em><sub>2</sub><em><sub>,σ</sub></em><sub>2</sub>).

To further motivate this setting, consider simulating <em><sub>n </sub></em>= 10<em><sub>,</sub></em>000 heights from the population of both males and females. Assume that males are distributed normal with mean <em><sub>µ</sub></em><sub>1 </sub>= 70[in] and standard deviation <em>σ</em><sub>1 </sub>= 3[in] and females are distributed normal with mean <em><sub>µ</sub></em><sub>2 </sub>= 64[in] and standard deviation <em><sub>σ</sub></em><sub>2 </sub>= 2<em><sub>.</sub></em>5[in]. Also assume that each distribution contributes equal mass, i.e., set the mixture parameter to <em><sub>δ </sub></em>= <em><sub>.</sub></em>5. The distribution of males is governed by

<em>f</em><em>,</em>

and the distribution of females is governed by

<em>f</em><em>.</em>

Below shows the pdf of <em><sub>f</sub></em><sub>1</sub>(<em><sub>x</sub></em>;<em><sub>µ</sub></em><sub>1</sub><em><sub>,σ</sub></em><sub>1</sub>), <em><sub>f</sub></em><sub>2</sub>(<em><sub>x</sub></em>;<em><sub>µ</sub></em><sub>2</sub><em><sub>,σ</sub></em><sub>2</sub>) and the mixture <em><sub>f</sub></em>(<em><sub>x</sub></em>) all on the same plot.

<table width="632">

 <tbody>

  <tr>

   <td width="632">x &lt;- <strong>seq</strong>(45,90,by=.05) n.x &lt;- <strong>length</strong>(x) f_1 &lt;- <strong>dnorm</strong>(x,mean=70,sd=3) f_2 &lt;- <strong>dnorm</strong>(x,mean=64,sd=2.5) f &lt;- <strong>function</strong>(x) { <strong>return</strong>(.5<strong>*</strong><strong>dnorm</strong>(x,mean=70,sd=3) <strong>+ </strong>.5<strong>*</strong><strong>dnorm</strong>(x,mean=64,sd=2.5))}plot_df &lt;- <strong>data.frame</strong>(x=<strong>c</strong>(x,x,x),Density=<strong>c</strong>(f_1,f_2,<strong>f</strong>(x)),Distribution=<strong>c</strong>(<strong>rep</strong>(“Male”,n.x),<strong>rep</strong>(“Female”,n.x),<strong>rep</strong>(“Mixture”,n.x)))<strong>library</strong>(ggplot2) <strong>ggplot</strong>(data = plot_df) <strong>+ </strong><strong>geom_line</strong>(mapping = <strong>aes</strong>(x = x, y = Density,color=Distribution))<strong>+ </strong><strong>labs</strong>(title = “Mixture of Normals”)</td>

  </tr>

 </tbody>

</table>

<h2>Mixture of Normals</h2>

<h3>Part I: Simulating a Mixture of Normals</h3>

The first goal is to simulate from the mixture distribution

<em>δf</em><sub>1</sub>(<em>x</em>;<em>µ</em><sub>1</sub><em>,σ</em><sub>1</sub>) + (1 <em>− δ</em>)<em>f</em><sub>2</sub>(<em>x</em>;<em>µ</em><sub>2</sub><em>,σ</em><sub>2</sub>)<em>,</em>

where <em><sub>µ</sub></em><sub>1 </sub>= 70<em><sub>,σ</sub></em><sub>1 </sub>= 3<em><sub>,µ</sub></em><sub>2 </sub>= 64<em><sub>,σ</sub></em><sub>2 </sub>= 2<em><sub>.</sub></em>5<em><sub>,δ </sub></em>= <em><sub>.</sub></em>5. We use the accept-reject algorithm to accomplish this task.

First we must choose the “easy to simulate” distribution <em><sub>g</sub></em>(<em><sub>x</sub></em>). For this problem choose <em><sub>g</sub></em>(<em><sub>x</sub></em>) to be a Cauchy distribution centered at 66 with scale parameter 7.

<table width="632">

 <tbody>

  <tr>

   <td width="632">g &lt;- <strong>function</strong>(x) { s=7 l=66 <strong>return</strong>(1<strong>/</strong>(pi<strong>*</strong>s<strong>*</strong>(1<strong>+</strong>((x<strong>–</strong>l)<strong>/</strong>s)<strong>^</strong>2)))}</td>

  </tr>

 </tbody>

</table>

<h4>Perform the following tasks</h4>

1) Identify a <strong><sub>suitable </sub></strong>value of <em><sub>alpha </sub></em>such that your envelope function <em><sub>e</sub></em>(<em><sub>x</sub></em>) satisfies

<em>f</em>(<em>x</em>) <em>≤ e</em>(<em>x</em>) = <em>g</em>(<em>x</em>)<em>/α, </em>where 0 <em>&lt; α &lt; </em>1<em>.</em>

Note that you must choose <em><sub>α </sub></em>so that <em><sub>e</sub></em>(<em><sub>x</sub></em>) is close to <em><sub>f</sub></em>(<em><sub>x</sub></em>). There is not one unique solution to this problem. The below plot shows how <em><sub>α </sub></em>= <em><sub>.</sub></em>20 creates an envelope function that is too large. Validate your choice of <em><sub>alpha </sub></em>with with a graphic similar to below.

<table width="632">

 <tbody>

  <tr>

   <td width="632"><em># Choose alpha </em>alpha &lt;- .20<em># Define envelope e(x) </em>e &lt;- <strong>function</strong>(x) { <strong>return</strong>(<strong>g</strong>(x)<strong>/</strong>alpha) }<em># Plot</em>x.vec &lt;- <strong>seq</strong>(30,100,by=.1) <strong>ggplot</strong>() <strong>+ </strong><strong>geom_line</strong>(mapping = <strong>aes</strong>(x = x.vec, y = <strong>f</strong>(x.vec)),col=”purple”)<strong>+ </strong><strong>geom_line</strong>(mapping = <strong>aes</strong>(x = x.vec, y = <strong>e</strong>(x.vec)),col=”green”)</td>

  </tr>

 </tbody>

</table>

<table width="632">

 <tbody>

  <tr>

   <td width="316">                                                40                                            60</td>

   <td width="316">                                        80                                           100x.vec</td>

  </tr>

  <tr>

   <td width="316"><em># Is g(x)&gt;f(x)?</em><strong>all</strong>(<strong>e</strong>(x.vec)<strong>&gt;</strong><strong>f</strong>(x.vec))</td>

   <td width="316"> </td>

  </tr>

 </tbody>

</table>

## [1] TRUE

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

2) Write a function named <strong><sub>r.norm.mix() </sub></strong>that simulates <strong><sub>n.samps </sub></strong>from the normal-mixture <em><sub>f</sub></em>(<em><sub>x</sub></em>). To accomplish this task you will wrap a function around the accept-reject algorithm from the lecture notes. Also include the acceptance rate, i.e., how many times did the algorithm accept a draw compared to the total number of trials performed. Your function should return a list of two elements: (i) the simulated vector mixture and (ii) the proportion of accepted cases. Run your function <strong><sub>r.norm.mix() </sub></strong>to simulate 10,000 cases and display the first 20 values. What’s the proportion of accepted cases? Compare this number to your chosen <em><sub>α </sub></em>and comment on the result. The code below should help you get started.

<h4>Solution</h4>

<table width="632">

 <tbody>

  <tr>

   <td width="632">r.norm.mix &lt;- <strong>function</strong>(n.samps) {<em>#n &lt;- 0 # counter for number samples accepted</em><em>#m &lt;- 0 # counter for number of trials</em><em>#samps                      &lt;- numeric(n.samps) # initialize the vector of output</em><em>#while (n &lt; n.samps) {</em><em># Fill in body ——-</em><em>#}</em><em>#return(list(x=samps,alpha.hat=n.samps/m))</em>}</td>

  </tr>

 </tbody>

</table>

<em>#r.norm.mix(n.samps=10000)$alpha.hat</em>

<em>#alpha</em>

3) Using <strong><sub>ggplot </sub></strong>or <strong><sub>base R</sub></strong>, construct a histogram of the simulated mixture distribution with the true mixture pdf <em><sub>f</sub></em>(<em><sub>x</sub></em>) overlayed on the plot.

<strong>Solution</strong>

<em>## Solution goes here ——-</em>

<h3>Part II: Bayesian Statistics and MCMC</h3>

Suppose that the experimenter collected 100 cases from the true mixture-normal distribution <em><sub>f</sub></em>(<em><sub>x</sub></em>). To solve problems (4) through (8) we analyze one realized sample from our function <strong><sub>r.norm.mix()</sub></strong>. In practice this dataset would be collected and not simulated. Uncomment the below code to simulate our dataset <strong><sub>x</sub></strong>. If you failed to solve Part I, then read in the csv file <strong><sub>mixture_data.csv </sub></strong>posted on Canvas.

<h4>Solution</h4>

<em># Simulate data</em>

<em>#set.seed(1983)</em>

<em>#x &lt;- r.norm.mix(n.samps=100)$x</em>

<em>#head(x)</em>

<em>#hist(x,breaks=20,xlab=”X”,main=””)</em>

<em># Or read data </em>x &lt;- <strong>read.csv</strong>(“mixture_data.csv”)<strong>$</strong>x <strong>head</strong>(x)

## [1] 71.66666 63.91096 67.06554 65.49516 70.34363 65.69982

<strong>hist</strong>(x,breaks=20,xlab=”X”,main=””)

X

Further, suppose that we know the true heights and standard deviations of the two normal distributions but the mixture parameter <em><sub>δ </sub></em>is unknown. In this case, we know <em><sub>µ</sub></em><sub>1 </sub>= 70, <em><sub>σ</sub></em><sub>1 </sub>= 3, <em><sub>µ</sub></em><sub>2 </sub>= 64, <em><sub>σ</sub></em><sub>2 </sub>= 2<em><sub>.</sub></em>5. The goal of this exercise is to utilize <strong>maximum likelihood </strong>and <strong>MCMC Bayesian techniques </strong>to estimate mixture parameter <em><sub>δ</sub></em>.

<h3>Maximum likelihood Estimator of Mixture Parameter</h3>

4) Set up the likelihood function <em><sub>L</sub></em>(<em><sub>δ</sub>|</em><em><sub>x</sub></em><sub>1</sub><em><sub>,…,x</sub></em><sub>100</sub>) and define it as <strong><sub>mix.like()</sub></strong>. The function should have two inputs including the parameter <strong><sub>delta </sub></strong>and data vector <strong><sub>x</sub></strong>. Evaluate the likelihood at the parameter values <strong><sub>delta=.2</sub></strong>, <strong><sub>delta=.4, </sub></strong>and <strong><sub>delta=.6</sub></strong>. Note that all three evaluations will be very small numbers. Which delta (<em><sub>δ </sub></em>= <em><sub>.</sub></em>2<em><sub>,.</sub></em>4<em><sub>,.</sub></em>6) is the most likely to have generated the dataset <strong><sub>x</sub></strong>?

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

5) Compute the maximum likelihood estimator of mixture parameter <em><sub>δ</sub></em>. To accomplish this task, apply your likelihood function <strong>mix.like() </strong>across the vector <strong>seq(.1,.99,by=.001)</strong>. The solution to this exercise is given below.

<em># delta &lt;- seq(.1,.99,by=.001)</em>

<em># MLE.values &lt;- sapply(delta,mix.like,x=x)</em>

<em># delta.MLE &lt;- delta[which.max(MLE.values)]</em>

<em># plot(delta,MLE.values,ylab=”Likelihood”,type=”l”)</em>

<em># abline(v=delta.MLE,col=”blue”)</em>

<em># text(x=delta.MLE+.08,y=mix.like(delta=.45,x=x),paste(delta.MLE),col=”blue”)</em>

<h3>MCMC</h3>

6) Run the Metropolis-Hastings algorithm to estimate mixture parameter <em><sub>δ</sub></em>. In this exercise you will assume a Beta(<em><sub>α </sub></em>= 10<em><sub>,β </sub></em>= 10) prior distribution on mixture parameter <em><sub>δ</sub></em>. Some notes follow:

<ul>

 <li>Run 20000 iterations. I.e., simulate 20000 draws of <em><sub>δ</sub></em><sup>(<em>t</em>)</sup></li>

 <li>Proposal distribution Beta(<em><sub>α </sub></em>= 10<em><sub>,β </sub></em>= 10) • Independence chain with Metropolis-Hastings ratio:</li>

</ul>

<em>R</em>(<em>δ</em>(<em>t</em>)<em>,δ</em><em>∗</em>) = <em>LL</em>((<em>δδ</em>(<em>∗t</em>)<em>||xx</em>11<em>,…,x,…,x</em>100100))

Display the first 20 simulated cases of <em><sub>δ</sub></em><sup>(<em>t</em>)</sup>.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

7) Construct a lineplot of the simulated Markov chain from exercise (6). The vertical axis is the simulated chain <em><sub>δ</sub></em><sup>(<em>t</em>) </sup>and the horizontal axis is the number of iterations.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

8) Plot the empirical autocorrelation function of your simulated chain <em><sub>δ</sub></em><sup>(<em>t</em>)</sup>. I.e., run the function <strong><sub>acf()</sub></strong>. A quick decay of the chain’s autocorrelations indicate good mixing properties.

<h4>Solution</h4>

<em>#acf(delta_t_vec,main=”ACF: Prior Beta(10,10)”)</em>

9) Compute the empirical Bayes estimate <em><sub>δ</sub></em>ˆ<em><sub>B </sub></em>of the simulated posterior distribution <em><sub>π</sub></em>(<em><sub>δ</sub>|</em><em><sub>x</sub></em><sub>1</sub><em><sub>,…,xn</sub></em>). To solve this problem, simply compute the sample mean of your simulated chain <em><sub>δ</sub></em><sup>(<em>t</em>) </sup>after discarding a 20% burn-in.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

10) Construct a histogram of the simulated posterior <em><sub>π</sub></em>(<em><sub>δ</sub>|</em><em><sub>x</sub></em><sub>1</sub><em><sub>,…,xn</sub></em>) after discarding a 20% burn-in.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

11) Run the Metropolis-Hastings algorithm to estimate the mixture parameter <em><sub>δ </sub></em>using a Beta(<em><sub>α </sub></em>= 15<em><sub>,β </sub></em>= 2) prior distribution on mixture parameter <em><sub>δ</sub></em>. Repeat exercises 6 though 10 using the updated prior.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

<strong>lineplot:</strong>

Construct a lineplot of the simulated Markov chain from exercise (6). The vertical axis is the simulated chain <em>δ</em><sup>(<em>t</em>) </sup>and the horizontal axis is the number of iterations.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

<strong>ACF:</strong>

Plot the empirical autocorrelation function of your simulated chain <em><sub>δ</sub></em><sup>(<em>t</em>)</sup>. I.e., run the function <strong><sub>acf()</sub></strong>. A slow decay of the chain’s autocorrelations indicate poor mixing properties.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

<strong>Bayes estimate:</strong>

Compute the empirical Bayes estimate <em><sub>δ</sub></em>ˆ<em><sub>B </sub></em>of the simulated posterior distribution <em><sub>π</sub></em>(<em><sub>δ</sub>|</em><em><sub>x</sub></em><sub>1</sub><em><sub>,…,xn</sub></em>). To solve this problem, simply compute the sample mean of your simulated chain <em><sub>δ</sub></em><sup>(<em>t</em>) </sup>after discarding a 20% burn-in. Your answer should be close to the MLE.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>

<strong>Posterior: </strong>Construct a histogram of the simulated posterior <em><sub>π</sub></em>(<em><sub>δ</sub>|</em><em><sub>x</sub></em><sub>1</sub><em><sub>,…,xn</sub></em>) after discarding a 20% burn-in.

<h4>Solution</h4>

<em>## Solution goes here ——-</em>