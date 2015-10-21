#### The Markov Switching (MS) model declares a relationship between a variable $y_{t}$ and a non observable variable $s_{t}$. The non observable variable can have one of $N$ different states. The model MS assumes that, once the past $m$ values of the non observable variable are known, the dependent variable follows the density $f(y_{t}|\mathcal{Y}_{t},s_{t}=j_{t},s_{t-1}=j_{t-1}...s_{t-m}=j_{t-m}\theta)$, where $\mathcal{Y}_{t-1}=\{y_{t-1},...y_{1}\}$. The states of the variable $s_{t}$ are determined by a discrete first order Markov Chain.

#### The model assumes the next probability distribution for $s_{t}$:

$$p(s_{t}=j|s_{t-1}=i,...s_{0}=i_{0},\mathbf{p})=p(s_{t}=j|s_{t-1}=i)=p_{ij}$$

#### In this program I implement the results found by @Hamilton90 and @Kim94 to estimate a family of models of the MS type using the Expected Maximization algorithm. The MS models has a wide range of specificatios, the most general specification that this program is able to estimate is 

$$\mathbf{y}_{t}=\mu_{s_{t}}+\phi_{s_{t}}(\mathbf{y}_{t-1}-\mu_{s_{t-1}})+...\phi_{s_{t}}(\mathbf{y}_{t-m}-\mu_{s_{t-m}})+\phi_{x,s_{t}}\mathbf{x}_{t}+e_{t},\ \ e_{t}\sim N(0,\Sigma_{s_{t}})\label{General}$$

#### other models can be estimated as a variation of ([General]). For example one can estimate the same equation as ([General]) but without switching in the varianze covarianze matrix $\Sigma$.Following @krolzig2013markov the model can be describe as MS<span>\*</span> where <span>\*</span> is the combination of the characters:

<span>M</span>
:   when the model has switching in the mean

<span>I</span>
:   when the model has switching in the intercept term

<span>A</span>
:   when the model has switching in the autoregressive terms

<span>X</span>
:   when the model has switching in the exogen variables

<span>H</span>
:   when the model has switching in the variance-covarianze matrix


