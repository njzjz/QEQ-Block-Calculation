# QEq Fragmentation Calculation study.
**Note: Fragmentation does not accelerate the QEq calculation.**

After running [calculateerror.py](calculateerror.py), we found QEq fragmentation is accurate.

However, we found fragmentation does not accelerate the QEq calculation. This may be because the time complexity of QEQ is **O(N)** since QEq uses a CG algorithm, so the fragmentation cannot accelerating its calculation.

**Author**: Jinzhe Zeng
**Email**: jzzeng@stu.ecnu.edu.cn
