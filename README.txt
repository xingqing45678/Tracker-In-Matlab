创建日期：2017-6-9
修改日期：2017-6-9

作者：qw
E-mail：1406820937@qq.com

内容：整合了跟踪算法的matlab代码，增加了原视频文件的兼容格式，可以做到代码性能间的对比

GitHub 操作流程：
1.  在本地编辑的工程目录文件夹下右键建立git init 初始化
2.  git status 查看当前目录下文件
3.  复制或者自己写.gitignore文件（touch .gitignore 命令创建该文件）
4.  git add -A 将文件目录下所有符合规定的文件加入（加入后可用git status查看加入成功与否）
5.  提交当前工作空间的修改内容，执行 git commit -m "提交信息" 将文件提交到repository里，提交信息用英文的双引号括起来。
6.  这时运行 git log 就可以看到提交的记录了
7.  在网页端创建好了远程仓库repository，先到Github上复制远程仓库的SSH地址。
8.  运行 git remote add origin 你复制的地址，如：git remote add origin git@github.com:xingqing45678/HashTracker.git 将本地库与远程库关联。
9.  执行 git push -u origin master 将本地仓库上传至Github的仓库并进行关联。
10.  以后想在commit后同步到Github上，只要直接执行 git push 就行了，可以在Github上看到修改。

每次修改本地文件后，在文件目录邮件建立github，然后执行4、5、10步凑即可
