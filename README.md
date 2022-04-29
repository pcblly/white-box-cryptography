# white-box-cryptography

**Author：CZX**

**Email：[chenzx1009@foxmail.com](mailto:chenzx1009@foxmail.com)**

> 《SM4分组密码白盒研究与实现》源代码

## Ⅰ. 实验环境

对于SM4分组密码白盒研究与实现的实验环境如下。

| **软件** | **版本**                                |
| -------- | --------------------------------------- |
| OS       | Ubuntu 20.04LTS                         |
| Kernel   | 5.13.0-30-generic                       |
| CPU      | Intel(R) Core(TM) i7-8550U CPU @1.80GHZ |
| G++      |                                         |
| NTL      | 11.5.1                                  |

> [NTL 是一个高性能、可移植的 C++ 库，为任意长度的整数提供数据结构和算法；用于整数和有限域上的向量、矩阵和多项式；以及任意精度的浮点运算。](https://libntl.org/doc/tour-intro.html)

NTL安装见[官方教程](https://libntl.org/doc/tour-unix.html)。

## Ⅱ. 运行方法

下载源码并解压源代码，进入`Project/src`目录。

**1.白盒复现**

```shell
#编译
g++ -g -O2 -std=c++11 -pthread -march=native main.cpp ../lib/*.cpp -o main -I ../include -I /usr/include -lntl -lgmp -lm
#运行
bash ./main
```

**2 攻击**

```shell
#编译
g++ -g -O2 -std=c++11 -pthread -march=native Attack_main.cpp ../lib/*.cpp -o attack -I ../include -I /usr/include -lntl -lgmp -lm
#运行
bash ./attack
```


