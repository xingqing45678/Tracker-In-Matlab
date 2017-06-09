This code was developed by Yuanhao Gong during his PhD at MOSAIC Group. 

Please cite Yuanhao's PhD thesis if you use this code in your work. Thank you!

=============================================================
@phdthesis{gong:phd,
  title={Spectrally regularized surfaces},
  author={Gong, Yuanhao},
  year={2015},
  school={ETH Zurich, Nr. 22616},
  note={http://dx.doi.org/10.3929/ethz-a-010438292}}
=============================================================

FAQ:

1) Why dual mesh (DM) structure is needed?
There are two reasons. First, these four sets guarantee the convergence. Second, 
we can use the updated neighbors for current position. Therefore, it is more computational efficient.

2) What is the difference between these three filters?
In general, GC filter is better in preserving details, compared with the other two. And
TV filter is better in removing noise as well as details. MC filter is between these two.

These three filters are correspond to three types of variational models. User should decide
which prior is to be assumed about the ground truth. 

3) What is the difference between split and nosplit scheme?
In general, splitting the image into four sets and looping on them is computational faster.
However, in some cases like deconvolution, we need to merge the four sets after every iteration.
So, it is better do nosplit scheme.

These two lead to exactly the same result. The split code is just more cache friendly.

此代码是由袁浩巩博士在叶组中发展起来的。
请举出Yuanhao的博士论文，如果你在工作中使用此代码。谢谢您!
=============================================================
“博士论文{龚：博士，
标题= {谱正则曲面}，
作者= {宫}，Yuanhao，
年份= { 2015 }，
苏黎世联邦理工大学学校= { 22616 }，NR，
注= { HTTP：/ / DX。DOI。org / 10.3929 / ethz-a-010438292 } }
=============================================================
常见问题解答：
1）为什么需要双网（干）结构？
有两个原因。首先，这四套保证收敛。二，
我们可以使用更新的邻居来进行当前位置。因此，它是更有效的计算。
2）这3个过滤器之间的区别是什么？
在一般情况下，GC过滤器是更好的保存细节。和TV过滤器是更好的消除噪音以及细节。MC滤波器介于这两个之间。
这三个过滤器是对应于三种类型的变分模型。要确定用户
哪一优先要根据实际情况。
3）分裂和不分裂方案之间的区别是什么？
在一般情况下，将图像分为四组，循环对他们是计算速度更快。
然而，在某些情况下，如反褶积，我们需要合并后，每一次迭代的四组。
所以，最好是nosplit方案。
这两个导致完全相同的结果。拆分代码只是更多的缓存友好。