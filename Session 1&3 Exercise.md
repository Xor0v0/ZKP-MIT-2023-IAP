# Session 1&3 Exercise

## 0x01 Session 1

##### Exercise 1

Question: Currently, you can only select adjacent pairs of nodes to check. Would the proof still be zero knowledge if you could pick arbitrary pairs of nodes to check?

答：协议不再是零知识的。原因：如果验证者可以选择任意一对节点，则会破坏协议的零知识性。具体而言，缺少了相邻节点颜色不同的约束，验证者总是可以从不相邻节点中获取知识，并且恢复出该图的完整三色染色解决方案。

##### Exercise 2

Question: The equation currently being used for confidence is `1-(1/E)^n`, where `E` is the number of edges in the graph, and `n` is the number of trials run. Is this the correct equation? Why is there no prior?

答：正确的。所谓prior，即先验，先验概率是根据先前的经验所得到的概率。由于在协议中每次挑选相邻节点来验证，事件（即挑选一条边）是相互独立的，所以并不存在先验。

##### Optional - ZKP for DLOG

Implement a non-interactive ZKP for discrete log in code! To do this, you’ll need to read and understand the first section of [this handout](https://people.eecs.berkeley.edu/~jfc/cs174/lecs/lec24/lec24.pdf), as well as the [Fiat-Shamir heuristic](https://en.wikipedia.org/wiki/Fiat–Shamir_heuristic).

其中涉及到了很多数学知识，比如mod运算、欧拉定理、DLOG等等，请参考本项目的学习笔记。

代码取自大佬的[github](https://github.com/dajuguan/MIT-IAP-2023/blob/main/Session_1_code.py)

当质数p取大一点之后，效率明显下降。具体而言，当p取16bit时，则会引发RuntimeError.wwod

##### zkmessage.xyz

Q1: Explain why you need to generate and save a “secret” value.

A1: Because I need a private key to generate a public key and then sign my message.

Q2: Write out a plain-English explanation of what statement is being proven in ZK.

A2: The author of this post is in the list but you don't know who is the real author.

Q3: Log into the same zkmessage account, from a different browser or computer. Explain why zkmessage can’t just use a simple “username/password” system like most social apps.

A3: Because "username/password" system need a centralized server to store user information，in that case, server can modify or even delete out post. 

## 0x02 Session 3