# DNA-Fountain illustrated(Adapted from dna-fountain)


## 总览

`demo.py`演示了将`test.txt`编码之后再恢复成原文件。

`encode.py`用于将原文件编码成一组DNA序列，`decode.py`用于解码恢复成原文件。

### 参数详解：

* 1 byte = 4 nt        因为一个碱基能编码两个比特，一个字节是8个比特。
* DNA序列长度是100nt， DNA序列 = 4 byte seed + 16 byte payload + 5 byte reedsolo
* payload size 和 chunk size 相同， payload由若干个chunk异或得到， 原文件按照chunk size分割成chunks，若构不成完整的chunk,则补零。

---

## 用法

> 先决条件

1. 安装 python3
2. 安装 numpy,reedsolo,tqdm。一种安装方式：`$pip3 install numpy reedsolo tqdm`
