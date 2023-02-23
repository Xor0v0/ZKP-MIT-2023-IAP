# 零知识证明｜part 1

通过零知识证明系列笔记，不仅要学到ZKP的原理及其应用，其中蕴含的思想和技术也可以散发到其他领域，比如可验证计算VC、可验证数据库VDB。希望有所收益！

共学网址：[zkp-co-learn](https://github.com/Antalpha-Labs/zkp-co-learn/discussions/2)

参考文章：

>  数学参考：[MIT数学讲义](https://mit6875.github.io/HANDOUTS/numbertheory.pdf)，[整数模N循环群](https://chenliang.org/2021/03/04/multiplicative-group-of-integers-modulo-n/)
>
> 认知ZKP参考：[zkp-intro](https://github.com/sec-bit/learning-zkp/tree/master/zkp-intro)

在part 1，主要做三件事：

- 数学基础
- 密码学基础
- 认知ZKP（Zero Knowledge Proof）

## 0x01 数学基础

#### 1. 初等数论

- 质数

	令$\pi(N)$作为小雨N的所有质数的数目，有：$\pi(N)=\Omega(N/\log N)$。因此，我们随机挑选一个小于N的正整数，有$1/\log N$的概率获得一个质数。

	如何高效判断一个数是否是素数？素性检验：AKS: M. Agrawal, N. Kayal, and N. Saxena. Primes is in p. *Annals of mathematics*, pages 781–793,2004.

	如何生成随机质数？如果我要生成`n-bit`的质数，只需要选择一个小于$2^n$随机数，判断其是否是质数，经过预期n轮循环，你将获得一个`n-bit`随机质数。

	个人而言，直接调python库：

	```python
	from crypto.Util.number import getPrime, isPrime
	```

	但是这两个函数内部实现没看，应该是安全的吧。

- 最大公因数｜gcd

	`gcd(a, b)`表示整除正整数a, b的最大正整数。一般采用欧几里得算法：

	```python
	def gcd(a, b) {
	  if a < b :
	  	a, b = b, a
	  while a % b != 0 :
	  	a, b =b, a % b
	  return b
	}
	```

- 扩展欧几里得算法｜ext_gcd

	扩展欧几里得算法主要是求两个整数x, y，满足`ax + by = gcd(a, b)`：

	```python
	def ext_gcd(a, b):
	  if a < b :
	    a, b = b, a
	   if b == 0:
	    	return 1, 0
	   else:
	      x, y = ext_gcd(b, a%b)
	      x, y = y, x - (a//b) * y
	      return x, y
	```

​		简单推导一下：当b==0时，显然x和y分别是1和0。当b!=0时，有：$gcd(a, b)=gcd(b, a\%b)$

​		对于左边，$gcd(a, b)=ax+by$，

​		对于右边，$gcd(b, a\%b)=bx_1+(a\%b)*y_1=bx_1+(a-(a/b)*b)*y_1=ay_1+b(x_1-(a/b)*y_1)$,

​		左右对应有：$x=y_1, y=(x_1-(a/b)*y_1)$。推毕。

​		该算法主要用于求乘法逆元。

- 模运算

	如果$ab=1(\mod N)$，则称b是a的模N乘法逆元，其存在条件是$gcd(a,N)=1$。

	运用ext_gcd算法即可求得乘法逆元：

	```python
	def get_inverse(a, N):
	  if gcd(a, N) == 1:
	    x, y = ext_gcd(a, N)
	    return x
	  print("No inverse!")
	  return 0
	```

	模幂运算：

	直接先幂后模是很低效计算方式，目前主要使用repeated squaring algorithm（时间复杂度O(n^3)）

	```python
	# 递归快速幂
	def qpow(x, exp):
	    if exp == 0:
	        return 1
	 	if exp % 2 == 0:
	        half = qpow(x, exp // 2)
	    	return half * half
	  	else:
	    	return x * qpow(x, exp - 1)
	```

- 中国剩余定理｜CRT

	CRT用于解决线性一元同余方程组问题。（所有的模数两两互质）
	$$
	\left\{\begin{align}
		x &= a_1 \mod n_1 \\
		x &= a_2 \mod n_2 \\
		&... \\
		x &= a_i \mod n_i
	\end{align}\right.
	$$
	

	用人话说一下原理：要求满足要求的`x`值，那么很朴素的想法是，在模$n_i$的时候，$x_i\equiv a_i\mod n_i$，此时模除$n_i$的数都等于0，即$x_i\equiv 0\mod n_j$。最后所有的情况$x_i$都相加即可。

	```python
	# 有物不知其数，三三数之剩二，五五数之剩三，七七数之剩二。问物几何？
	from gmpy2 import invert
	import numpy as np
	
	# a: array, n: array, k: int
	def CRT(a, n, k):
	    nn = np.prod(n)
	    ans = 0
	    for i in range(k):
	        m = nn // n[i]
	        m_ = invert(int(m), int(n[i]))
	        ans += a[i] * m * m_
	    return ans % nn
	
	if __name__ == '__main__':
	    a, n = np.array([2, 3, 2]), np.array([3, 5, 7])
	    print(CRT(a, n, 3))
	```

- 欧拉函数｜$\varphi(N)$

	$\varphi(N)$表示在小于N的正整数中与N互质的个数。

	由于任意大于1的正整数都可以分解成若干质数的乘积，不妨设$N=\prod_{i=1}^lp_i^{\alpha_i}$，则$\varphi(N)=\prod_{i=1}^lp_i^{\alpha_i-1}(p_i-1)$。特别地，如果$N=p^\alpha$，则$\varphi(N)=p^{\alpha-1}(p-1))$。更为特殊地，如果$N=p$，则$\varphi(N)=p-1$。

- 欧拉定理

	对于任意两个正整数a和n，若两者互质，则有：$a^{\varphi(n)}\equiv 1\mod N$。特殊情况下，若n为质数，记作p，有$a^{p-1}\equiv 1\mod n$，这个式子也被称为费马小定理。

- 二次剩余｜Quadratic Residues

	

#### 2. 群论｜Groups

##### 关于群的基本概念

- 一般群：群是由一个元素集合以及定义在集合上的运算构成。它满足封闭性（closure）、单位元（Identity）、逆元（inverse）、结合律（Associativity）。
- 单位元｜$\mathcal{I}$：群里任何元素与单位元作运算得到的结果是仍然是该元素。
- 逆元：群里的任何元素与其逆元做运算得到的结果是单位元。
- 交换群或者阿贝尔群｜Abelian Group：除了遵循一般群的法则，还要求满足交换律（Commmutativity）。

- 群阶和群元素的阶：群阶就是群中元素的个数；群元素的阶是指该元素变为单位元所需要执行运算的次数。
- 拉格朗日定理｜Lagrange's Theorem：群内任意元素的阶整除群的阶。
- 生成元｜generator：一个群的生成元就是阶为群阶的元素，记为g，群可表示为$\mathbb{G}=\{g,g^2,\dots,g^{|\mathbb{G}|}=\mathcal{I}\}$
- 循环群｜Cyclic group：具有生成元的群叫做循环群。群$\mathbb{Z}_N$就是一个循环群。**注意，循环群的生成元不唯一。**
- 定理2：如果一个群的阶为素数，则这个群是循环群，并且群中除了单位元都是该群的生成元。
- 离散对数｜DLOG：令$\mathbb{G}$是一个循环群，g为其生成器，那么对于群的任一元素我们可以记作$h=g^x$。定义$x=d\log_g(h)$，则称x是h以g为底的离散对数。**在密码学中，我们需要寻找计算h很容易，但是计算离散对数很难的这样的群。**群$\mathbb{Z}_N$就是一个计算离散对数容易的群，这是因为$\mathbb{Z}_N$的生成器就是与N互质的那些数，如果你知道了$x\cdot g=h\mod N$，那么你就可以利用扩展欧几里得算法算出g的乘法逆元，进而算出x，所以$\mathbb{Z}_N$不可用。**群$\mathbb{Z}_N^*$被猜测是计算dlog困难的群**。

##### 整数模N乘法群$\mathbb{Z}_N^*$

在同余理论中，模N的**互质同余类**组成一个乘法群，称为整数模 N 乘法群。从定义我们可以知道，群$\mathbb{Z}_N^*$的阶就是$\varphi(N)$

定理3：$\mathbb{Z}_N^*$是一个乘法运算和求逆运算都能在多项式时间内完成的群。

##### 整数模P乘法群$\mathbb{Z}_P^*$（P为素数）

群$\mathbb{Z}_P^*$的阶的为$\varphi(P)=P-1$.

定理4: 如果P是素数，则$\mathbb{Z}_P^*$是一个循环群。

- 如何计算群$\mathbb{Z}_P^*$的生成元数目？

	定理5: 群$\mathbb{Z}_P^*$的生成元数目是$\varphi(P-1)$.

- 如何高效的判断给定的元素是否是生成元？

	令给定元素为x，则可计算$x^{|\mathbb{Z}|}\equiv 1\mod P$，若成立则是生成元，否则不是。

- 如何对群$\mathbb{Z}_P^*$的生成元进行随机抽样？

#### 3. 椭圆曲线｜Eliptic Curve

椭圆曲线的数学形式：$y^2=x^3+Ax+B$

需要知道的是：椭圆曲线与有限域上的群有对应关系

#### 4. 多项式｜polynomail



## 0x02 密码学基础

#### 1. 安全性定义

- 香农给出了密码系统的完美保密性（perfect secret）的定义，即攻击者无法从密文中获取明文的任何信息。我们称之为完美安全，或者说无条件安全。
- 完美安全的要求太高，人们对其进行了弱化relax：攻击者在多项式时间内无法从密文中获取明文的任何信息。
- 进一步弱化：提出了计算安全、整体安全等安全性定义。

#### 2. 安全假设

正如RSA密码体制建立在大整数分解难题上（但量子计算机可以攻破，Shor算法），所有密码体制均是建立在数学难题上，本小节介绍ZKP协议所依赖的常见的数学难题。**请注意，这些数学难题只是目前的理论水平上证明甚至是猜测是困难的，但不意味着永远是困难的。由于量子计算机的出现，现有很多密码体制均会面临被攻破的风险，为此人们提出了后量子密码学用以利用新的数学难题建立密码体制，此处不做扩展**

##### DLOG难题

##### CDH难题

##### DDH难题

...

#### 3. 安全性证明

任何一个密码体制都需要对其进行安全性证明。在安全性证明过程中，首先要提出安全模型（所谓安全模型，其实就是一种假设，假设越强，其安全性越弱）。目前一般使用**一般群模型**和**随机预言模型**，二者并无强弱之分，使用场景不同罢了。

可证明安全的一般步骤

> 1. 首先判断密码体制$\mathcal{S}_1$所依赖的数学难题$\mathcal{P}$或者依赖的是其他的密码体制$\mathcal{S}_2$。
> 2. 其次，密码体制$\mathcal{S}_1$的攻击者$\mathcal{A}$构造所依赖数学难题的求解器$\mathcal{B}$或者其他的密码体制$\mathcal{S}_2$的攻击者$\mathcal{B}$。要求$\mathcal{B}$和$\mathcal{S}_1$的表现行为对于$\mathcal{A}$是不可区分的。

#### 4. 安全目标+攻击模型

不同的密码体制拥有不同的安全目标，研究人员讨论在不同的攻击模型下密码体制是否能够满足安全目标。

常见的安全目标有：

- 单向性：由密文无法推断出明文。
- 语义安全：由密文无法获知明文的任何信息。（**语义安全就是不可区分安全**）
- 不可区分性：敌手选择两个明文，密码体制随机选择加密其中一个，敌手无法从密文中判断加密的是哪个明文。
- 非延展性：敌手无法从密文中推断得到一个关于明文的有意义变换的密文。

常见的攻击模型有：

- 选择明文攻击CPA：允许攻击者访问预言机
- 非适应性选择密文攻击CCA1：攻击者在得到挑战密文之前，允许访问预言机
- 适应性选择密文攻击CCA2：攻击者在得到挑战密文之后，仍然允许访问预言机

## 0x03 认知ZKP

Zero Knowledge Proof，是由Goldwasser、Micali和Rackoff三位提出的一种运行在prover和verifier之间的两方密码协议，其运行机制一般是**挑战-响应机制（challenge-response）**，即下文介绍的sigma协议。首先弄清楚：

> 1. 什么是知识？
>
>    参考文章中总结的非常到位，醍醐灌顶。**知识是能够帮助人们提高计算能力的东西，与公共所知的东西有关。**
>
> 2. 什么是证明？
>
>    证明就是prover利用证据witness来使得verfier相信它拥有知识的过程。在计算机领域，证明的过程与计算密不可分。**需要注意：证明与验证往往是不对称的，证明往往需要消费大量的算力，而验证则比较简单，由此衍生出了可验证计算、云存储的完整性检验、可验证数据库等技术**。人们利用云具有强大算力和海量存储资源的特性，将计算服务、文件存储服务、数据库存储服务等外包给云，由其生成证明，数据拥有者则只需要验证结果即可。证明的过程可分为交互式的和非交互式的。

零知识证明ZKP技术，是prover在不向verifier泄漏知识的情况下，使其向证明其拥有知识。其具有三个特性：

- 完备性Completeness：只有当prover拥有知识时，它才能总是正确回答verifier的质询
- 可靠性Soundness：如果prover不拥有知识，它总会被逮住（即以不可忽略的概率无法正确回复质询）
- 零知识性Zero Knowledge：Prover的响应不会泄漏知识。

> 3. 那么如何证明协议具有零知识性呢？
>
>    参考文献引入了模拟的概念。其实我感觉跟安全性证明中可证明安全的随机预言模型有点像。。。当现实世界中拥有知识的Alice能够像理想世界中不拥有知识的Zlice一样，准确地回答Bob的每一次质询时，那么Bob无法区分此二者，即可以说协议具有零知识性。
>
>    在上述表述中，需要注意的是：现实世界的Alice无法具有「超能力」，否则这个协议将不具有零知识性。
>
> 4. 完备性和可靠性之间的对称性：完备性保证了诚实的Alice一定能成功，可靠性保证了恶意的Alice一定会失败。
>
> 5. 关于协议可靠性的证明，参考插图非常精美直观的参考文章。

##### Schnorr协议（Σ-协议）

协议的建立依赖于离散对数数学难题：给定有限域上的整数`r`，在循环群中找到对应的点`r*G`是简单的，但反过来计算却是很难的。

协议流程如下（图片来源于参考文章）：

![img](https://github.com/sec-bit/learning-zkp/raw/master/zkp-intro/3/img/schnorr.png)

整个流程清晰简单，需要注意的点：

- Alice选取随机数r的原因：保护知识（也就是私钥）无法被攻击者抽取。
- 随机数的生成需要具有密码学安全强度。否则将导致Bob将可以恢复出私钥，协议的可靠性就丧失掉了。

关于ZK更多术语、概念，我将在Part 2详细阐述、加以区分。