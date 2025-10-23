"""
Python Learn - Day 1 within 100 Days
Author: guisongwu
Date  : 2025-10-19
"""
# Guido van Rossum（吉多·范罗苏姆）
# Netherlands (荷兰)




"""
Python Learn - Day 2 within 100 Days - First Python Program
Author: guisongwu
Date  : 2025-10-19
"""
# This is a single line comment. Multiline comment use """  """
# print('hello world')
# print('hello world\n')
# print("hello world")




"""
Python Learn - Day 3 within 100 Days - Variables in Python
Author: guisongwu
Date  : 2025-10-19
"""
# # int
# print(0b100)  # 二进制整数(binary) 
# print(0o100)  # 八进制整数(octal)
# print(100)    # 十进制整数
# print(0x100)  # 十六进制整数(hexadecimal)

# # float
# print(123.456)    # 数学写法
# print(1.23456e2)  # 科学计数法

# # str
# print('hello')
# print("world")

# # bool
# print(True)
# print(False)

# # basic calc
# a = 45
# b = 12
# print(a, b)   # 45 12
# print(a + b)  # 57
# print(a - b)  # 33
# print(a * b)  # 540
# print(a / b)  # 3.75

# # types
# print(type(100))  # <class 'int'>
# print(type(123.45))  # <class 'float'>
# print(type('h'))  # <class 'str'> There is not class <char> in python
# print(type('hello'))  # <class 'str'>
# print(type(True))  # <class 'bool'>

# # type invert
# a = 100 # int
# b = 123.45 # float
# c = '123' # str
# d = '100' # str
# e = '123.45' # str
# f = 'hello, world' # str
# g = True # bool
# print(float(a))         # int类型的100转成float，输出100.0
# print(int(b))           # float类型的123.45转成int，输出123
# print(int(c))           # str类型的'123'转成int，输出123
# print(int(c, base=16))  # str类型的'123'按十六进制转成int，输出291
# print(int(d, base=2))   # str类型的'100'按二进制转成int，输出4
# print(float(e))         # str类型的'123.45'转成float，输出123.45
# print(bool(f))          # str类型的'hello, world'转成bool，输出True. str类型转成bool类型时，只要字符串有内容，不是''或""，对应的布尔值都是True
# print(int(g))           # bool类型的True转成int，输出1
# print(chr(a))           # int类型的100转成str，输出'd'
# print(type(chr(a))) 
# print(ord('d'))         # str类型的'd'转成int，输出100




"""
Python Learn - Day 4 within 100 Days - Operators
Author: guisongwu
Date  : 2025-10-23
"""
# [] 索引
# [:] 切片
# ** 幂
# & 按位与
# ~ 按位取反
# % 取余
# // 取商

# print(321 // 12)    # 整除运算，输出26
# print(321 % 12)     # 求模运算，输出9
# print(321 ** 12)    # 求幂运算，输出1196906950228928915420617322241

# 赋值运算符
# SyntaxError: invalid syntax
# print((a = 10))
# := 海象运算符, 运算符右侧的值也是整个表达式的值
# print((a := 10))  # 10
# print(a)          # 10

# 逻辑运算符 and, or, not
# (False) and (), 则右边括号里的表达式不会被执行(短路)
# (True) or (), 右边括号里的表达式不会被执行(短路)
# 比较运算符的优先级高于赋值运算符，所以下面的flag0 = 1 == 1先做1 == 1产生布尔值True，再将这个值赋值给变量flag0。print函数可以输出多个值，多个值之间可以用,进行分隔，输出的内容默认以空格分开。
# flag0 = 1 == 1 # True
# flag1 = 3 > 2  # True
# flag2 = 2 < 1  # False
# flag3 = flag1 and flag2 # False
# flag4 = flag1 or flag2  # True
# flag5 = not flag0 # False
# print('flag0 =', flag0)     # flag0 = True
# print('flag1 =', flag1)     # flag1 = True
# print('flag2 =', flag2)     # flag2 = False
# print('flag3 =', flag3)     # flag3 = False
# print('flag4 =', flag4)     # flag4 = True
# print('flag5 =', flag5)     # flag5 = False
# print(flag1 and not flag2)  # True
# print(1 > 2 or 2 == 3)      # False

# 将华氏温度转换为摄氏温度
# f = float(input('请输入华氏温度: '))
# c = (f - 32) / 1.8
# print('%.1f华氏度 = %.1f摄氏度' % (f, c))
# print(f'{f:.1f}华氏度 = {c:.1f}摄氏度') # 字符串前面的f表示这个字符串是需要格式化处理的，其中的{f:.1f}和{c:.1f}可以先看成是{f}和{c}，表示输出时会用变量f和变量c的值替换掉这两个占位符，后面的:.1f表示这是一个浮点数，小数点后保留1位有效数字。

# r = float(input('请输入半径：'))
# area = 3.1415926 * r * r
# perimeter = 2 * 3.1415926 * r
# print('area = %.2f' % area)
# print(f'perimeter = {perimeter:.2f}')

# import math
# r = float(input('请输入半径：'))
# area = math.pi * r * r
# perimeter = 2 * math.pi * r
# print('area = %.2f' % area)
# print(f'perimeter = {perimeter:.2f}') # perimeter = 18.85
# print(f'{perimeter = }') # perimeter = 18.84955592153876


# year = int(input('请输入年份: '))
# is_leap = year % 4 == 0 and year % 100 != 0 or year % 400 == 0
# print(f'{is_leap = }')




"""
Python Learn - Day 5 within 100 Days - Judgment and branch
Author: guisongwu
Date  : 2025-10-23
"""










