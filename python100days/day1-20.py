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
Python Learn - Day 5 within 100 Days - Judgment and Banch
Author: guisongwu
Date  : 2025-10-23
"""
# BMI
# height = float(input('身高(cm)：'))
# weight = float(input('体重(kg)：'))
# bmi = weight / (height / 100) ** 2
# print(f'{bmi = :.1f}')
# print(f'height = {height:.2f} weight = {weight:.2f} bmi = {bmi:.2f}')
# if 18.5 <= bmi < 24:
#     print('你的身材很棒！')
# if bmi < 18.5:
#     print('你的体重过轻!')
# elif bmi < 24:
#     print('你的身材很棒!')
# elif bmi < 27:
#     print('你的体重过重!')
# else:
#     print('你的体重太重!')

# status_code = int(input('响应状态码: '))
# if status_code == 400:
#     description = 'Bad Request'
# elif status_code == 401:
#     description = 'Unauthorized'
# elif status_code == 403:
#     description = 'Forbidden'
# elif status_code == 404:
#     description = 'Not Found'
# elif status_code == 405:
#     description = 'Method Not Allowed'
# elif status_code == 418:
#     description = 'I am a teapot'
# elif status_code == 429:
#     description = 'Too many requests'
# else:
#     description = 'Unknown status Code'
# print('状态码描述:', description)

# status_code = int(input('响应状态码: '))
# match status_code:
#     case 400: description = 'Bad Request'
#     case 401: description = 'Unauthorized'
#     case 403: description = 'Forbidden'
#     case 404: description = 'Not Found'
#     case 405: description = 'Method Not Allowed'
#     case 418: description = 'I am a teapot'
#     case 429: description = 'Too many requests'
#     case _: description = 'Unknown Status Code' # 带有_的case语句在代码中起到通配符的作用，如果前面的分支都没有匹配上，代码就会来到case _。case _的是可选的，并非每种分支结构都要给出通配符选项。如果分支中出现了case _，它只能放在分支结构的最后面，如果它的后面还有其他的分支，那么这些分支将是不可达的。
# print('状态码描述:', description)

# status_code = int(input('响应状态码: '))
# match status_code:
#     case 400 | 405: description = 'Invalid Request'
#     case 401 | 403 | 404: description = 'Not Allowed'
#     case 418: description = 'I am a teapot'
#     case 429: description = 'Too many requests'
#     case _: description = 'Unknown Status Code'
# print('状态码描述:', description)

# 3x-5 x>1
# x+2 -1<=x<=1
# 5x+3 x<-1
# x = float(input('x: '))
# if x < -1:
#     y = 5*x+3
# elif -1 <= x <= 1:
#     y = x + 2
# else:
#     y = 3*x-5
# print(f"{y = }")
# print(f"y = {y:.2f}")

# score = float(input("请输入成绩: "))
# if score >= 90:
#     grade = 'A'
# elif score >= 80:
#     grade = 'B'
# elif score >= 70:
#     grade = 'C'
# elif score >= 60:
#     grade = 'D'
# else:
#     grade = 'E'
# print(f"{grade = }")

# a = float(input("a: "))
# b = float(input("b: "))
# c = float(input("c: "))
# if a + b > c and a + c > b and b + c > a:
#     perimeter = a + b + c
#     print(f"{perimeter = }")
#     s = perimeter / 2
#     area = (s * (s - a) * (s - b) * (s - c)) ** 0.5
#     print(f"{area = }")
# else:
#     print("This is not a triangle!")
    



"""
Python Learn - Day 6 within 100 Days - Iteration
Author: guisongwu
Date  : 2025-10-23
"""
# import time
# print("hello")
# time.sleep(1)
# print("world")

# 如果明确知道循环执行的次数，我们推荐使用 for-in 循环
# import time
# for i in range(5):
# for _ in range(5):
#     print("hello world")
#     time.sleep(1)

# print(type(range(100))) # class 'range'
# range(101)：可以用来产生0到100范围的整数，需要注意的是取不到101。
# range(1, 101)：可以用来产生1到100范围的整数，相当于是左闭右开的设定，即[1, 101)。
# range(1, 101, 2)：可以用来产生1到100的奇数，其中2是步长（跨度），即每次递增的值，101取不到。
# range(100, 0, -2)：可以用来产生100到1的偶数，其中-2是步长（跨度），即每次递减的值，0取不到。

# total = 0
# for i in range(1, 101):
#     total += i
# print(f"{total = }")

# total = 0
# for i in range(1, 101):
#     if i % 2 == 0:
#         total += i
# print(f"{total = }")

# total = 0
# for i in range(2, 101, 2):
#     total += i
# print("total = %d" % total)
# print(f"{total = }")
# print(f"total = {total}")
# print("sum = ", sum(range(2, 101, 2)))

# 要构造循环结构但是又不能确定循环重复的次数，我们推荐使用 while 循环
# total = 0
# i = 1
# while i <= 100:
#     total += i
#     i += 1
# print(total)

# break and continue
# total = 0
# i = 2
# while True:
#     total += i
#     i += 2
#     if i > 100:
#         break
# print(total) 

# total = 0
# for i in range(1, 101):
#     if i % 2 != 0:
#         continue
#     total += i
# print(total)

# for i in range(1, 10):
#     for j in range(1, i + 1):
#         print(f'{i}×{j}={i * j}', end='\t') # end 参数控制每次 print() 输出后结尾要打印的内容, 默认是"\n"
#     print()

# judge if a num is prime or not
# num = int(input("give a positive num: "))
# end = int(num ** 0.5)
# is_prime = True
# for i in range(1, end+1):
#     if num % i == 0:
#         is_prime = False
#         break
# if is_prime:
#     print(f"{num} is a prime")
# else:
#     print(f"{num} is not a prime")

# 最大公因数
# a = int(input("a = "))
# b = int(input("b = "))
# for i in range(a, 0, -1):
#     if a % i == 0 and b % i == 0:
#         print(f"最大公因数{i}")
#         break

# 辗转相除
# x = int(input('x = '))
# y = int(input('y = '))
# while y % x != 0:
#     x, y = y % x, x
# print(f'最大公约数: {x}')')

# 猜数字
# import random
# answer = random.randrange(0, 100)
# counter = 0
# while True:
#     guess = int(input("your guess: "))
#     counter += 1
#     if guess < answer:
#         print("Greater!")
#     elif guess == answer:
#         print("Bingo!")
#         break
#     else:
#         print("Smaller!")
# print(f"Success after {counter} times!")




"""
Python Learn - Day 7 within 100 Days - Practice with Branching and Loop Structures
Author: guisongwu
Date  : 2025-10-24
"""
# 100 以内的素数
# for i in range(2, 100):
#     is_prime = True
#     for j in range(2, int(i**0.5)+1):
#         if i % j == 0:
#             is_prime = False
#             break # break只能终止它所在的那个循环，这一点在使用嵌套循环结构时需要引起注意
#     if is_prime:
#         print(i)

# Fibonacci Sequence
# a, b = 0, 1
# for i in range(20):
#     a, b = b, a+b
#     print(a)
# “同时解包赋值”（tuple unpacking assignment）
# 右边先整体求值，生成一个元组：(b, a + b), 使用的是赋值前的 a 和 b 的旧值
# 再一次性将结果解包并赋给左边的变量：
# a = 旧的b
# b = 旧的a + 旧的b

# narcissistic number (水仙花数)
# for i in range(100, 1000):
#     a = i % 10
#     b = ((i - a) // 10) % 10
#     c = i // 100
#     if i == (a**3 + b**3 + c**3):
#         print(i)

# invert a num
# num = int(input("num = "))
# inverted_num = 0
# while num > 0:
#     inverted_num = inverted_num * 10 + num % 10
#     num = num // 10
#     print(inverted_num) 
#     print(num) 

# 公鸡 5 元一只，母鸡 3 元一只，小鸡 1 元三只，用 100 块钱买一百只鸡，问公鸡a、母鸡b、小鸡c各有多少只？
# for a in range(0, 21):
#     for b in range(0, 34):
#         for c in range(0, 100, 3):
#             if a + b + c == 100 and 5*a + 3*b + 1./3.*c == 100:
#                 print(a, b, c)
#                 break

# CRAPS赌博游戏. 说明：CRAPS又称花旗骰，是美国拉斯维加斯非常受欢迎的一种的桌上赌博游戏。该游戏使用两粒骰子，玩家通过摇两粒骰子获得点数进行游戏。简化后的规则是：玩家第一次摇骰子如果摇出了 7 点或 11 点，玩家胜；玩家第一次如果摇出 2 点、3 点或 12 点，庄家胜；玩家如果摇出其他点数则游戏继续，玩家重新摇骰子，如果玩家摇出了 7 点，庄家胜；如果玩家摇出了第一次摇的点数，玩家胜；其他点数玩家继续摇骰子，直到分出胜负。为了增加代码的趣味性，我们设定游戏开始时玩家有 1000 元的赌注，每局游戏开始之前，玩家先下注，如果玩家获胜就可以获得对应下注金额的奖励，如果庄家获胜，玩家就会输掉自己下注的金额。游戏结束的条件是玩家破产（输光所有的赌注）。
# import random
# money = 1000
# while money > 0:
#     a = int(input("请下注："))
#     first_toss = random.randrange(1, 7) + random.randrange(1, 7)
#     print(f"first toss {first_toss}")
#     if first_toss == 7 or first_toss == 11:
#         money += a
#         print(f"You won {a}")
#     elif first_toss == 2 or first_toss == 3 or first_toss == 12:
#         money -= a
#         print(f"You lose {a}")
#     else:
#         while True:
#             second_toss = random.randrange(1, 7) + random.randrange(1, 7)
#             print(f"second toss {second_toss}")
#             if second_toss == 7:
#                 money -= a
#                 print(f"You lose {a}")
#                 break
#             elif second_toss == first_toss:
#                 money += a
#                 print(f"You won {a}")
#                 break
#     print(f"Now your money is {money}")
    



"""
Python Learn - Day 8 within 100 Days - List 1
Author: guisongwu
Date  : 2025-10-24
"""
# 容器型变量
# items1 = [35, 12, 99, 68, 55, 35, 87]
# items2 = ['Python', 'Java', 'Go', 'Kotlin']
# items3 = [100, 12.3, 'Python', True]
# print(items1)  # [35, 12, 99, 68, 55, 35, 87]
# print(type(items1))  # <class 'list'>
# print(items2)  # ['Python', 'Java', 'Go', 'Kotlin']
# print(items3)  # [100, 12.3, 'Python', True]

# 通过 Python 内置的list函数将其他序列变成列表。准确的说，list并不是一个普通的函数，它是创建列表对象的构造器
# print(type(range(1, 10))) # <class 'range'>
# items4 = list(range(1, 10))
# print(type('hello')) # <class 'str'>
# items5 = list('hello')
# print(items4)  # [1, 2, 3, 4, 5, 6, 7, 8, 9]
# print(items5)  # ['h', 'e', 'l', 'l', 'o']

# 可以使用+运算符实现两个列表的拼接，拼接运算会将两个列表中的元素连接起来放到一个列表中
# 可以使用*运算符实现列表的重复运算，*运算符会将列表元素重复指定的次数
# 可以使用in或not in运算符判断一个元素在不在列表中
# items5 = [35, 12, 99, 45, 66]
# items6 = [45, 58, 29]
# items7 = ['Python', 'Java', 'JavaScript']
# print(items5 + items6)  # [35, 12, 99, 45, 66, 45, 58, 29]
# print(items6 + items7)  # [45, 58, 29, 'Python', 'Java', 'JavaScript']
# items5 += items6
# print(items5)  # [35, 12, 99, 45, 66, 45, 58, 29]

# print(items6 * 3)  # [45, 58, 29, 45, 58, 29, 45, 58, 29]
# print(items7 * 2)  # ['Python', 'Java', 'JavaScript', 'Python', 'Java', 'JavaScript']

# print(29 in items6)  # True
# print(99 in items6)  # False
# print('C++' not in items7)     # True
# print('Python' not in items7)  # False

# 当我们想操作列表中的某个元素时，可以使用[]运算符，通过在[]中指定元素的位置来访问该元素，这种运算称为索引运算。需要说明的是，[]的元素位置可以是0到N - 1的整数，也可以是-1到-N的整数，分别称为正向索引和反向索引，其中N代表列表元素的个数。对于正向索引，[0]可以访问列表中的第一个元素，[N - 1]可以访问最后一个元素；对于反向索引，[-1]可以访问列表中的最后一个元素，[-N]可以访问第一个元素，代码如下所示。
# 避免出现索引越界的情况，对于下面的items8，如果我们访问items8[5]或items8[-6]，就会引发IndexError错误，导致程序崩溃，对应的错误信息是：list index out of range
# items8 = ['apple', 'waxberry', 'pitaya', 'peach', 'watermelon']
# print(items8[0])   # apple
# print(items8[2])   # pitaya
# print(items8[4])   # watermelon
# items8[2] = 'durian'
# print(items8)      # ['apple', 'waxberry', 'durian', 'peach', 'watermelon']
# print(items8[-5])  # 'apple'
# print(items8[-4])  # 'waxberry'
# print(items8[-1])  # watermelon
# items8[-4] = 'strawberry'
# print(items8)      # ['apple', 'strawberry', 'durian', 'peach', 'watermelon']
# # 希望一次性访问列表中的多个元素，我们可以使用切片运算。切片运算是形如[start:end:stride]的运算符，其中start代表访问列表元素的起始位置，end代表访问列表元素的终止位置（终止位置的元素无法访问）
# print(items8[1:3:1])     # ['strawberry', 'durian']
# print(items8[0:3:1])     # ['apple', 'strawberry', 'durian']
# print(items8[0:5:2])     # ['apple', 'durian', 'watermelon']
# print(items8[-4:-2:1])   # ['strawberry', 'durian']
# print(items8[-2:-6:-1])  # ['peach', 'durian', 'strawberry', 'apple']
# # 如果start值等于0，那么在使用切片运算符时可以将其省略；如果end值等于N，N代表列表元素的个数，那么在使用切片运算符时可以将其省略；如果stride值等于1，那么在使用切片运算符时也可以将其省略。所以，下面的代码跟上面的代码作用完全相同。
# print(items8[1:3])     # ['strawberry', 'durian']
# print(items8[:3:1])    # ['apple', 'strawberry', 'durian']
# print(items8[:3:])    # ['apple', 'strawberry', 'durian']
# print(items8[::2])     # ['apple', 'durian', 'watermelon']
# print(items8[-4:-2])   # ['strawberry', 'durian']
# print(items8[-2::-1])  # ['peach', 'durian', 'strawberry', 'apple']
# # 切片修改列表
# items8[1:3] = ['x', 'o']
# print(items8)  # ['apple', 'x', 'o', 'peach', 'watermelon']

# 两个列表还可以做关系运算，我们可以比较两个列表是否相等，也可以给两个列表比大小
# nums1 = [1, 2, 3, 4]
# nums2 = list(range(1, 5))
# nums3 = [3, 2, 1]
# print(nums1 == nums2)  # True
# print(nums1 != nums2)  # False
# print(nums1 <= nums3)  # True
# print(nums2 >= nums3)  # False

# 如果想逐个取出列表中的元素，可以使用for-in循环
# method 1
# languages = ["python", "c", "c++", "matlab"]
# print(len(languages))
# for i in range(len(languages)):
#     print(languages[i])
# # method 2
# for language in languages:
#     print(language)

# toss a dice
# import random
# counters = [0] * 6
# # 模拟掷色子记录每种点数出现的次数
# for _ in range(100000):
#     face = random.randrange(1, 7)
#     counters[face - 1] += 1
#     # 输出每种点数出现的次数
# for face in range(1, 7):
#     print(f'{face}点出现了{counters[face - 1]}次')



"""
Python Learn - Day 9 within 100 Days - List 2
Author: guisongwu
Date  : 2025-10-24
"""
# Methods of list
# items.append(item)
# items.insert(index, item)

# items.remove(item)
# items.pop()
# items.pop(index)
# items.clear()

# items.index(item, index)
# items.count(item)

# items.sort()
# items.reverse()


# 列表的方法
# 列表是一种可变容器，可变容器指的是我们可以向容器中添加元素、可以从容器移除元素，也可以修改现有容器中的元素。我们可以使用列表的append方法向列表中追加元素，使用insert方法向列表中插入元素。追加指的是将元素添加到列表的末尾，而插入则是在指定的位置添加新元素
# languages = ['Python', 'Java', 'C++']
# languages.append('JavaScript')
# print(languages)  # ['Python', 'Java', 'C++', 'JavaScript']
# languages.insert(1, 'SQL')
# print(languages)  # ['Python', 'SQL', 'Java', 'C++', 'JavaScript']
# languages.insert(2, 'Matlab')
# print(languages)  # ['Python', 'SQL', 'Matlab', 'Java', 'C++', 'JavaScript']

# 我们可以用列表的remove方法从列表中删除指定元素，需要注意的是，如果要删除的元素并不在列表中，会引发ValueError错误导致程序崩溃，所以建议大家在删除元素时，先用之前讲过的成员运算做一个判断。我们还可以使用pop方法从列表中删除元素，pop方法默认删除列表中的最后一个元素，当然也可以给一个位置，删除指定位置的元素。在使用pop方法删除元素时，如果索引的值超出了范围，会引发IndexError异常，导致程序崩溃。除此之外，列表还有一个clear方法，可以清空列表中的元素
# languages = ['Python', 'SQL', 'Java', 'C++', 'JavaScript']
# if 'Java' in languages:
#     languages.remove('Java')
# if 'Swift' in languages:
#     languages.remove('Swift')
# print(languages)  # ['Python', 'SQL', C++', 'JavaScript']
# languages.pop()
# print(languages)  # ['Python', 'SQL', C++']
# temp = languages.pop(1)
# print(temp)       # SQL
# print(type(temp))       # <class 'str'>
# languages.append(temp)
# print(languages)  # ['Python', C++', 'SQL']
# languages.clear()
# print(languages)  # []

# a = 1
# del a
# print(a)

# languages = ["Python", 'Matlab', 'Python']
# print(languages)
# if 'Python' in languages:
#     languages.remove('Python') # remove the first 'Python'
# print(languages)

# 元素位置和频次. 列表的index方法可以查找某个元素在列表中的索引位置，如果找不到指定的元素，index方法会引发ValueError错误；列表的count方法可以统计一个元素在列表中出现的次数
# items = ['Python', 'Java', 'Java', 'C++', 'Kotlin', 'Python']
# print(items.index('Python'))     # 0
# print(items.index('Python', 1))  # 5  从索引位置1开始查找'Python'
# print(items.count('Python'))     # 2
# print(items.count('Kotlin'))     # 1
# print(items.count('Swfit'))      # 0
# # print(items.index('Java', 3))    # ValueError: 'Java' is not in list  从索引位置3开始查找'Java'
# print(items.index('C++', 0))

# items = ['Python', 'Java', 'C++', 'Kotlin', 'Swift']
# items.sort()
# print(items)  # ['C++', 'Java', 'Kotlin', 'Python', 'Swift']
# items.reverse()
# print(items)  # ['Swift', 'Python', 'Kotlin', 'Java', 'C++']

# 列表生成式. 在 Python 中，列表还可以通过一种特殊的字面量语法来创建, 强烈建议用生成式语法来创建列表
# 场景一：创建一个取值范围在1到99且能被3或者5整除的数字构成的列表。
# method 1
# items = []
# for i in range(1, 100):
#     if i % 3 == 0 or i % 5 == 0:
#         items.append(i)
# print(items)
# method 2
# items = [i for i in range(1, 100) if i % 3 == 0 or i % 5 == 0]
# print(items)

# 场景二：有一个整数列表nums1，创建一个新的列表nums2，nums2中的元素是nums1中对应元素的平方。
# method 1
# nums1 = [1, 2, 3, 4, 5, 12, 14]
# nums2 = []
# for num in nums1:
#     nums2.append(num ** 2)
# print(nums2)
# method 2
# nums1 = [1, 2, 3, 4, 5, 12, 14]
# nums2 = [num ** 2 for num in nums1]
# print(nums2)

# 场景三： 有一个整数列表nums1，创建一个新的列表nums2，将nums1中大于50的元素放到nums2中。
# method 1
# nums1 = [35, 12, 97, 64, 55]
# nums2 = []
# for num in nums1:
#     if num > 50:
#         nums2.append(num)
# print(nums2)
# method 2
# nums1 = [35, 12, 97, 64, 55]
# nums2 = [num for num in nums1 if num > 50]
# print(nums2)


# 嵌套列表
# scores = [[95, 83, 92], [80, 75, 82], [92, 97, 90], [80, 78, 69], [65, 66, 89]]
# print(scores[0])
# print(scores[0][1])

# import random
# scores = [[random.randrange(60, 101) for _ in range(3)] for _ in range(5)]
# print(scores)

# import random
# red_balls = list(range(1, 34))
# selected_balls = []
# for _ in range(6):
#     index = random.randrange(len(red_balls))
#     # 将选中的球从红色球列表中移除并添加到选中列表
#     selected_balls.append(red_balls.pop(index))
# selected_balls.sort()
# for ball in selected_balls:
#     print(f'\033[031m{ball:0>2d}\033[0m', end=' ')
# blue_ball = random.randrange(1, 17)
# print(f'\033[034m{blue_ball:0>2d}\033[0m')
# 上面代码中print(f'\033[0m...\033[0m')是为了控制输出内容的颜色，红色球输出成红色，蓝色球输出成蓝色。其中省略号代表我们要输出的内容，\033[0m是一个控制码，表示关闭所有属性，也就是说之前的控制码将会失效，你也可以将其简单的理解为一个定界符，m前面的0表示控制台的显示方式为默认值，0可以省略，1表示高亮，5表示闪烁，7表示反显等。在0和m的中间，我们可以写上代表颜色的数字，比如30代表黑色，31代表红色，32代表绿色，33代表黄色，34代表蓝色等。]]])

# import random
# red_balls = [i for i in range(1, 34)]
# blue_balls = list(range(1, 17))
# # 从红色球列表中随机抽出6个红色球（无放回抽样）
# selected_balls = random.sample(red_balls, 6)
# selected_balls.sort()
# for ball in selected_balls:
#     print(f'\033[031m{ball:0>2d}\033[0m', end=' ')
# # 从蓝色球列表中随机抽出1个蓝色球
# blue_ball = random.choice(blue_balls)
# print(f'\033[034m{blue_ball:0>2d}\033[0m')

# import random
# from rich.console import Console
# from rich.table import Table 
# # table 是用 Rich 库（一个漂亮的终端输出库）创建的表格对象
# # 创建控制台
# console = Console()
# n = int(input('生成几注号码: '))
# red_balls = list(range(1, 34))
# blue_balls = [i for i in range(1, 17)]
# # 创建表格并添加表头
# table = Table(title='balls', show_header=True)
# for col_name in ['序号', '红球', '蓝球']:
#     table.add_column(col_name, justify='center')
# for i in range(n):
#     selected_balls = random.sample(red_balls, 6)
#     selected_balls.sort()
#     blue_ball = random.choice(blue_balls)
#     # 向表格中添加行（序号，红色球，蓝色球）
#     table.add_row(str(i + 1), 
#                   f'[red]{" ".join([f"{ball:0>2d}" for ball in selected_balls])}[/red]',
#                   f'[blue]{blue_ball:0>2d}[/blue]')
# # 通过控制台输出表格
# console.print(table)

# rich, table, console
# from rich.table import Table
# from rich.console import Console
# table = Table(title="Scores", show_header=True)
# console = Console()
# table.add_column("Name", justify="center")
# table.add_column("Score", justify = "center")
# table.add_row("xiaohong", "90")
# table.add_row("xiaoming", "80")
# console.print(table)


# a = [1, 3, 5]
# print(" ".join([str(i) for i in a]))



"""
Python Learn - Day 10 within 100 Days - Tuple
Author: guisongwu
Date  : 2025-10-25
"""
# 元组和列表的不同之处在于，元组是不可变类型，这就意味着元组类型的变量一旦定义，其中的元素不能再添加或删除，而且元素的值也不能修改。如果试图修改元组中的元素，将引发TypeError错误，导致程序崩溃.
# t1 = (35, 12, 98)
# t2 = ('jason', 24, True, 'Beijing')

# print(type(t1))  # <class 'tuple'>
# print(type(t2))  # <class 'tuple'>
# print(len(t1))  # 3
# print(len(t2))  # 4
# print(t1[0])    # 35
# print(t1[2])    # 98
# print(t2[-1])   # Beijing

# # 切片运算 [start:end:stride]
# print(t2[:2])   # ('jason', 24)
# print(t2[::3])  # ('jason', 'Beijing')
# for elem in t1:
#     print(elem)

# # 成员运算
# print(12 in t1)         # True
# print(99 in t1)         # False
# print('Hao' not in t2)  # True

# # 拼接运算
# t3 = t1 + t2
# print(t3)  # (35, 12, 98, '骆昊', 43, True, '四川成都')
# # 比较运算
# print(t1 == t3)            # False
# print(t1 >= t3)            # False
# print(t1 <= (35, 11, 99))  # False


# 如果元组中只有一个元素，需要加上一个逗号，否则()就不是代表元组的字面量语法，而是改变运算优先级的圆括号
# a = ()
# print(type(a))  # <class 'tuple'>
# b = ('hello')
# print(type(b))  # <class 'str'>
# c = (100)
# print(type(c))  # <class 'int'>
# d = ('hello', )
# print(type(d))  # <class 'tuple'>
# e = (100, )
# print(type(e))  # <class 'tuple'>

# 打包和解包操作. 当我们把多个用逗号分隔的值赋给一个变量时，多个值会打包成一个元组类型；当我们把一个元组赋值给多个变量时，元组会解包成多个值然后分别赋给对应的变量
# a = 1, 10, 100
# print(type(a))  # <class 'tuple'>
# print(a)        # (1, 10, 100)
# # 解包操作
# i, j, k = a
# print(i, j, k)  # 1 10 100

# 在解包时，如果解包出来的元素个数和变量个数不对应，会引发ValueError异常，错误信息为：too many values to unpack（解包的值太多）或not enough values to unpack（解包的值不足）
# a = 1, 10, 100, 1000
# i, j, k = a             # ValueError: too many values to unpack (expected 3)
# i, j, k, l, m, n = a    # ValueError: not enough values to unpack (expected 6, got 4)
# 有一种解决变量个数少于元素的个数方法，就是使用星号表达式。通过星号表达式，我们可以让一个变量接收多个值，代码如下所示。需要注意两点：
# 首先，用星号表达式修饰的变量会变成一个列表，列表中有0个或多个元素；其次，在解包语法中，星号表达式只能出现一次。
# a = 1, 10, 100, 1000
# i, j, *k = a
# print(i, j, k)        # 1 10 [100, 1000]
# print(type(i))
# print(type(j))
# print(type(k))
# i, *j, k = a
# print(i, j, k)        # 1 [10, 100] 1000
# *i, j, k = a
# print(i, j, k)        # [1, 10] 100 1000
# *i, j = a
# print(i, j)           # [1, 10, 100] 1000
# i, *j = a
# print(i, j)           # 1 [10, 100, 1000]
# i, j, k, *l = a
# print(i, j, k, l)     # 1 10 100 [1000]
# i, j, k, l, *m = a
# print(i, j, k, l, m)  # 1 10 100 1000 []

# 需要说明一点，解包语法对所有的序列都成立，这就意味着我们之前讲的列表、range函数构造的范围序列甚至字符串都可以使用解包语法。
# 用星号修饰的变量会变成一个列表 list
# a, b, *c = range(1, 10)
# print(a, b, c)
# print(type(a)) # <class 'int'>
# print(type(b)) # <class 'int'>
# print(type(c)) # <class 'list'>
# a, b, c = [1, 10, 100]
# print(a, b, c)
# a, *b, c = 'hello'
# print(a, b, c)
# print(type(a))
# print(type(b))
# print(type(c))

# 元组是不可变类型，不可变类型更适合多线程环境，因为它降低了并发访问变量的同步化开销。通常不可变类型在创建时间上优于对应的可变类型。我们可以使用timeit模块的timeit函数来看看创建保存相同元素的元组和列表各自花费的时间，timeit函数的number参数表示代码执行的次数。下面的代码中，我们分别创建了保存1到9的整数的列表和元组，每个操作执行10000000次，统计运行时间
# import timeit
# print('%.3f 秒' % timeit.timeit('[1, 2, 3, 4, 5, 6, 7, 8, 9]', number=10000000))
# print('%.3f 秒' % timeit.timeit('(1, 2, 3, 4, 5, 6, 7, 8, 9)', number=10000000))
# print('%.3f seconds' % timeit.timeit('[1,2,3,4,5,6]', number = 10000000))
# print('%.3f seconds' % timeit.timeit('(1,2,3,4,5,6)', number = 10000000))

# list and tuple can be transformed
# infos = ('骆昊', 43, True, '四川成都')
# print(list(infos))  # ['骆昊', 43, True, '四川成都']
# frts = ['apple', 'banana', 'orange']
# print(tuple(frts))  # ('apple', 'banana', 'orange')

# 列表和元组都是容器型的数据类型，即一个变量可以保存多个数据，而且它们都是按一定顺序组织元素的有序容器。列表是可变数据类型，元组是不可变数据类型，所以列表可以添加元素、删除元素、清空元素、排序反转，但这些操作对元组来说是不成立的。列表和元组都可以支持拼接运算、成员运算、索引运算、切片运算等操作，后面我们要讲到的字符串类型也支持这些运算，因为字符串就是字符按一定顺序构成的序列，在这一点上三者并没有什么区别。我们推荐大家使用列表的生成式语法来创建列表，它不仅好用而且效率很高，是 Python 语言中非常有特色的语法。



"""
Python Learn - Day 11 within 100 Days - String
Author: guisongwu
Date  : 2025-10-26
"""
# s1 = 'hello, world!'
# s2 = "你好，世界！❤️"
# s3 = '''hello,
# wonderful
# world!'''
# print(s1)
# print(s2)
# print(s3)
# s1 = '\'hello, world!\''
# s2 = '\\hello, world!\\'
# print(s1)
# print(s2)

# Python 中有一种以r或R开头的字符串，这种字符串被称为原始字符串，意思是字符串中的每个字符都是它本来的含义，没有所谓的转义字符。例如，在字符串'hello\n'中，\n表示换行；而在r'hello\n'中，\n不再表示换行，就是字符\和字符n
# s1 = 'it is \time \to \read \now' # \r是回车符（carriage return）相当于让输出回到了行首
# s2 = r'\it \is \time \to \read \now'
# print(s1)
# print(s2)

# +, *, [], [:], in, not in...
# s1 = 'hello' + ', ' + 'world'
# print(s1)    # hello, world
# s2 = '!' * 3
# print(s2)    # !!!
# s1 += s2
# print(s1)    # hello, world!!!
# s1 *= 2
# print(s1)    # hello, world!!!hello, world!!!

# s1 = 'a whole new world'
# s2 = 'hello world'
# print(s1 == s2)             # False
# print(s1 < s2)              # True ord('a') < ord('h')
# print(s1 == 'hello world')  # False
# print(s2 == 'hello world')  # True
# print(s2 != 'Hello world')  # True
# s3 = '骆昊'
# print(ord('骆'))            # 39558
# print(ord('昊'))            # 26122
# s4 = '王大锤'
# print(ord('王'))            # 29579
# print(ord('大'))            # 22823
# print(ord('锤'))            # 38180
# print(s3 >= s4)             # True
# print(s3 != s4)             # True

# print(ord('a'))
# print(chr(97))

# 成员运算符 in / not in
# s1 = 'hello, world'
# s2 = 'goodbye, world'
# print('wo' in s1)      # True
# print('wo' not in s2)  # False
# print(s2 in s1)        # False

# list1 = list(range(10))
# tuple1 = tuple(range(10))
# str1 = 'hello, world'
# print(len(list1))                 # 10
# print(len(tuple1))                 # 10
# print(len(str1))                 # 12

# string 不是可变类型
# s = 'abc123456'
# n = len(s)
# print(s[0], s[-n])    # a a
# print(s[n-1], s[-1])  # 6 6
# print(s[2], s[-7])    # c c
# print(s[5], s[-4])    # 3 3
# print(s[2:5])         # c12
# print(s[-7:-4])       # c12
# print(s[2:])          # c123456
# print(s[:2])          # ab
# print(s[::2])         # ac246
# print(s[::-1])        # 654321cba
# print(type(s[0]))
# print(type(s[3]))

# s = 'hello'
# for i in range(len(s)-1):
#     print(s[i], end="")
# print(s[-1])
# for elem in s:
#     print(elem)

# .capitalize(), .title(), .upper(), .lower()
# s1 = 'hello, world!'
# # 字符串首字母大写
# print(s1.capitalize())  # Hello, world!
# # 字符串每个单词首字母大写
# print(s1.title())       # Hello, World!
# # 字符串变大写
# print(s1.upper())       # HELLO, WORLD!
# s2 = 'GOODBYE'
# # 字符串变小写
# print(s2.lower())       # goodbye
# print(s1)               # hello, world
# print(s2)               # GOODBYE
# s3 = '123456'
# print(s3.upper())       # HELLO, WORLD!
# 字符串是不可变类型，使用字符串的方法对字符串进行操作会产生新的字符串，但是原来变量的值并没有发生变化。所以上面的代码中，当我们最后检查s1和s2两个变量的值时，s1和s2 的值并没有发生变化。


# 查找操作. 如果想在一个字符串中从前向后查找有没有另外一个字符串，可以使用字符串的find或index方法。在使用find和index方法时还可以通过方法的参数来指定查找的范围，也就是查找不必从索引为0的位置开始。
# find() / rfind()：找不到时返回 -1, 
# index() / rindex()：找不到时抛出 ValueError
# s = 'hello, world!'
# print(s.find('or'))      # 8
# print(s.find('or', 9))   # -1 没找到
# print(s.find('of'))      # -1
# print(s.index('or'))     # 8
# print(s.index('or', 9))  # ValueError: substring not found
# find和index方法还有逆向查找（从后向前查找）的版本，分别是rfind和rindex
# s = 'hello world!'
# print(s.find('o'))       # 4
# print(s.rfind('o'))      # 7
# print(s.rindex('o'))     # 7
# print(s.rindex('o', 8))  # ValueError: substring not found

# s1 = 'hello, world!'
# print(s1.startswith('He'))   # False
# print(s1.startswith('hel'))  # True
# print(s1.endswith('!'))      # True
# s2 = 'abc123456'
# print(s2.isdigit())  # False
# print(s2.isalpha())  # False
# print(s2.isalnum())  # True

# 格式化, 在 Python 中，字符串类型可以通过center、ljust、rjust方法做居中、左对齐和右对齐的处理。如果要在字符串的左侧补零，也可以使用zfill方法。
# s = 'hello, world'
# print(s.center(20, '*'))  # ****hello, world****
# print(s.rjust(20))        #         hello, world
# print(s.ljust(20, '~'))   # hello, world~~~~~~~~
# print('33'.zfill(5))      # 00033
# print('-33'.zfill(5))     # -0033
# print('abc'.zfill(5))     # -0033

# a = 321
# b = 123
# print('%d * %d = %d' % (a, b, a * b))
# # 当然，我们也可以用字符串的format方法来完成字符串的格式，代码如下所示。
# print('{0} * {1} = {2}'.format(a, b, a * b))
# # 从 Python 3.6 开始，可以在字符串前加上f来格式化字符串，{变量名}是一个占位符，会被变量对应的值将其替换掉，代码如下所示。
# print(f'{a} * {b} = {a * b}')

# 变量值          占位符      格式化结果           说明
# 3.1415926       {:.2f}        '3.14'          保留小数点后两位
# 3.1415926       {:+.2f}       '+3.14'         带符号保留小数点后两位
# -1              {:+.2f}       '-1.00'         带符号保留小数点后两位
# 3.1415926       {:.0f}        '3'             不带小数
# 123             {:0>10d}      '0000000123'    左边补0，补够10位
# 123             {:x<10d}      '123xxxxxxx'    右边补x ，补够10位
# 123             {:>10d}       '       123'    左边补空格，补够10位
# 123             {:<10d}       '123       '    右边补空格，补够10位
# 123456789       {:,}          '123,456,789'   逗号分隔格式
# 0.123           {:.2%}        '12.30%'        百分比格式
# 123456789       {:.2e}        '1.23e+08'      科学计数法格式

# 字符串的strip方法可以帮我们获得将原字符串修剪掉左右两端指定字符之后的字符串，默认是修剪空格字符。这个方法非常有实用价值，可以用来将用户输入时不小心键入的头尾空格等去掉，strip方法还有lstrip和rstrip两个版本，相信从名字大家已经猜出来这两个方法是做什么用的。
# s1 = '   jackfrued@126.com  '
# print(s1.strip())      # jackfrued@126.com
# print(s1)
# s2 = '~你好，世界~'
# print(s2.lstrip('~'))  # 你好，世界~
# print(s2.rstrip('~'))  # ~你好，世界
# print(s2)  # ~你好，世界

# 替换操作. 如果希望用新的内容替换字符串中指定的内容，可以使用replace方法，代码如下所示。replace方法的第一个参数是被替换的内容，第二个参数是替换后的内容，还可以通过第三个参数指定替换的次数。
# s = 'hello, good world'
# print(s.replace('o', '@'))     # hell@, g@@d w@rld
# print(s)     # hell@, g@@d w@rld
# print(s.replace('o', '@', 1))  # hell@, good world

# 拆分与合并. 可以使用字符串的split方法将一个字符串拆分为多个字符串（放在一个列表中），也可以使用字符串的join方法将列表中的多个字符串连接成一个字符串，代码如下所示。
# s = 'I love you'
# words = s.split()
# print(words)            # ['I', 'love', 'you']
# print('~'.join(words))  # I~love~you

# 需要说明的是，split方法默认使用空格进行拆分，我们也可以指定其他的字符来拆分字符串，而且还可以指定最大拆分次数来控制拆分的效果，代码如下所示。
# s = 'I#love#you#so#much'
# words = s.split('#')
# print(words)  # ['I', 'love', 'you', 'so', 'much']
# words = s.split('#', 2)
# print(words)  # ['I', 'love', 'you#so#much']

# 编码和解码. Python 中除了字符串str类型外，还有一种表示二进制数据的字节串类型（bytes）。
# 所谓字节串，就是由零个或多个字节组成的有限序列。
# 通过字符串的encode方法，我们可以按照某种编码方式将字符串编码为字节串，我们也可以使用字节串的decode方法，将字节串解码为字符串，代码如下所示。
# a = '骆昊'
# b = a.encode('utf-8')
# c = a.encode('gbk')
# print(b)                  # b'\xe9\xaa\x86\xe6\x98\x8a'
# print(type(b))            # <class 'bytes'>
# print(c)                  # b'\xc2\xe6\xea\xbb'
# print(type(c))            # b'\xe9\xaa\x86\xe6\x98\x8a'
# print(b.decode('utf-8'))  # 骆昊
# print(c.decode('gbk'))    # 骆昊

# 注意，如果编码和解码的方式不一致，会导致乱码问题（无法再现原始的内容）或引发UnicodeDecodeError错误，导致程序崩溃。



"""
Python Learn - Day 12 within 100 Days - Set
Author: guisongwu
Date  : 2025-10-26
"""
# 无序性说明集合中的元素并不像列中的元素那样存在某种次序，可以通过索引运算就能访问任意元素，集合并不支持索引运算。另外，集合的互异性决定了集合中不能有重复元素，这一点也是集合区别于列表的地方，我们无法将重复的元素添加到一个集合中。集合类型必然是支持in和not in成员运算的，这样就可以确定一个元素是否属于集合，也就是上面所说的集合的确定性。集合的成员运算在性能上要优于列表的成员运算，这是集合的底层存储特性决定的
# 创建集合. 在 Python 中，创建集合可以使用{}字面量语法，{}中需要至少有一个元素，因为没有元素的{}并不是空集合而是一个空字典，字典类型我们会在下一节课中为大家介绍。当然，也可以使用 Python 内置函数set来创建一个集合，准确的说set并不是一个函数，而是创建集合对象的构造器，这个知识点会在后面讲解面向对象编程的地方为大家介绍。我们可以使用set函数创建一个空集合，也可以用它将其他序列转换成集合，例如：set('hello')会得到一个包含了4个字符的集合（重复的字符l只会在集合中出现一次）。除了这两种方式，还可以使用生成式语法来创建集合，就像我们之前用生成式语法创建列表那样。
# set1 = {1, 2, 3, 3, 3, 2}
# print(set1)
# set2 = {'banana', 'pitaya', 'apple', 'apple', 'banana', 'grape'}
# print(set2)
# set3 = set('hello')
# print(set3)
# set4 = set([1, 2, 2, 3, 3, 3, 2, 1])
# print(set4)
# set5 = {num for num in range(1, 20) if num % 3 == 0 or num % 7 == 0}
# print(set5)
# 需要提醒大家，集合中的元素必须是hashable类型，所谓hashable类型指的是能够计算出哈希码的数据类型，通常不可变类型都是hashable类型，如整数（int）、浮点小数（float）、布尔值（bool）、字符串（str）、元组（tuple）等。可变类型都不是hashable类型，因为可变类型无法计算出确定的哈希码，所以它们不能放到集合中。例如：我们不能将列表作为集合中的元素；同理，由于集合本身也是可变类型，所以集合也不能作为集合中的元素。我们可以创建出嵌套列表（列表的元素也是列表），但是我们不能创建出嵌套的集合，这一点在使用集合的时候一定要引起注意。

# 元素的遍历. 我们可以通过len函数来获得集合中有多少个元素，但是我们不能通过索引运算来遍历集合中的元素，因为集合元素并没有特定的顺序。当然，要实现对集合元素的遍历，我们仍然可以使用for-in循环，代码如下所示。
# set1 = {'Python', 'C++', 'Java', 'Kotlin', 'Swift'}
# for elem in set1:
#     print(elem)
# 提示：大家看看上面代码的运行结果，通过单词输出的顺序体会一下集合的无序性。每次打印顺序都不一样

# 集合的运算. Python 为集合类型提供了非常丰富的运算，主要包括：成员运算、交集运算、并集运算、差集运算、比较运算（相等性、子集、超集）等。
# 成员运算. 可以通过成员运算in和not in 检查元素是否在集合中，代码如下所示。
# set1 = {11, 12, 13, 14, 15}
# print(10 in set1)      # False 
# print(15 in set1)      # True
# set2 = {'Python', 'Java', 'C++', 'Swift'}
# print('Ruby' in set2)  # False
# print('Java' in set2)  # True

# 二元运算. 集合的二元运算主要指集合的交集、并集、差集、对称差等运算，这些运算可以通过运算符来实现，也可以通过集合类型的方法来实现，代码如下所示。
# set1 = {1, 2, 3, 4, 5, 6, 7}
# set2 = {2, 4, 6, 8, 10}
# # 交集
# print(set1 & set2)                      # {2, 4, 6}
# print(set1.intersection(set2))          # {2, 4, 6}
# # 并集
# print(set1 | set2)                      # {1, 2, 3, 4, 5, 6, 7, 8, 10}
# print(set1.union(set2))                 # {1, 2, 3, 4, 5, 6, 7, 8, 10}
# # 差集
# print(set1 - set2)                      # {1, 3, 5, 7}
# print(set1.difference(set2))            # {1, 3, 5, 7}
# # 对称差
# print(set1 ^ set2)                      # {1, 3, 5, 7, 8, 10}
# print(set1.symmetric_difference(set2))  # {1, 3, 5, 7, 8, 10}

# 通过上面的代码可以看出，对两个集合求交集，&运算符和intersection方法的作用是完全相同的，使用运算符的方式显然更直观且代码也更简短。需要说明的是，集合的二元运算还可以跟赋值运算一起构成复合赋值运算，例如：set1 |= set2相当于set1 = set1 | set2，跟|=作用相同的方法是update；set1 &= set2相当于set1 = set1 & set2，跟&=作用相同的方法是intersection_update，代码如下所示。
# set1 = {1, 3, 5, 7}
# set2 = {2, 4, 6}
# # set1 |= set2
# set1.update(set2)
# print(set1)  # {1, 2, 3, 4, 5, 6, 7}
# set3 = {3, 6, 9}
# # set1 &= set3
# set1.intersection_update(set3)
# print(set1)  # {3, 6}
# # set2 -= set1
# set2.difference_update(set1)
# print(set2)  # {2, 4}

# 比较运算. 两个集合可以用==和!=进行相等性判断，如果两个集合中的元素完全相同，那么==比较的结果就是True，否则就是False。如果集合A的任意一个元素都是集合B的元素，那么集合A称为集合B的子集，即对于 ∀ a ∈ A ，均有 a ∈ B ，则 A ⊆ B ，A是B的子集，反过来也可以称B是A的超集。如果A是B的子集且A不等于B，那么A就是B的真子集。Python 为集合类型提供了判断子集和超集的运算符，其实就是我们非常熟悉的<、<=、>、>=这些运算符。当然，我们也可以通过集合类型的方法issubset和issuperset来判断集合之间的关系，代码如下所示。
# set1 = {1, 3, 5}
# set2 = {1, 2, 3, 4, 5}
# set3 = {5, 4, 3, 2, 1}
# print(set1 < set2)   # True
# print(set1 <= set2)  # True
# print(set2 < set3)   # False
# print(set2 <= set3)  # True
# print(set2 > set1)   # True
# print(set2 == set3)  # True
# print(set1.issubset(set2))    # True
# print(set2.issuperset(set1))  # True
# 说明：上面的代码中，set1 < set2判断set1是不是set2的真子集，set1 <= set2判断set1是不是set2的子集，set2 > set1判断set2是不是set1的超集。当然，我们也可以通过set1.issubset(set2)判断set1是不是set2的子集；通过set2.issuperset(set1)判断set2是不是set1的超集。

# 集合的方法. 刚才我们说过，Python 中的集合是可变类型，我们可以通过集合的方法向集合添加元素或从集合中删除元素。
# set1 = {1, 10, 100}
# # 添加元素
# set1.add(1000)
# set1.add(10000)
# print(set1)  # {1, 100, 1000, 10, 10000}
# # 删除元素
# set1.discard(10)
# if 100 in set1:
#     set1.remove(100)
# a = set1.pop()
# print(f"{a:0>10d}")
# print(set1)  # {1, 1000, 10000}
# # 清空元素
# set1.clear()
# print(set1)  # set()
# 说明：删除元素的remove方法在元素不存在时会引发KeyError错误，所以上面的代码中我们先通过成员运算判断元素是否在集合中。集合类型还有一个pop方法可以从集合中随机删除一个元素，该方法在删除元素的同时会返回（获得）被删除的元素，而remove和discard方法仅仅是删除元素，不会返回（获得）被删除的元素。

# 集合类型还有一个名为isdisjoint的方法可以判断两个集合有没有相同的元素，如果没有相同元素，该方法返回True，否则该方法返回False，代码如下所示。
# set1 = {'Java', 'Python', 'C++', 'Kotlin'}
# set2 = {'Kotlin', 'Swift', 'Java', 'Dart'}
# set3 = {'HTML', 'CSS', 'JavaScript'}
# print(set1.isdisjoint(set2))  # False
# print(set1.isdisjoint(set3))  # True

# 不可变集合. Python 中还有一种不可变类型的集合，名字叫frozenset。set跟frozenset的区别就如同list跟tuple的区别，frozenset由于是不可变类型，能够计算出哈希码，因此它可以作为set中的元素。除了不能添加和删除元素，frozenset在其他方面跟set是一样的，下面的代码简单的展示了frozenset的用法。
# fset1 = frozenset({1, 3, 5, 7})
# fset1 = frozenset([1, 3, 5, 7])
# fset2 = frozenset(range(1, 6))
# print(fset1)          # frozenset({1, 3, 5, 7})
# print(fset2)          # frozenset({1, 2, 3, 4, 5})
# print(fset1 & fset2)  # frozenset({1, 3, 5})
# print(fset1 | fset2)  # frozenset({1, 2, 3, 4, 5, 7})
# print(fset1 - fset2)  # frozenset({7})
# print(fset1 < fset2)  # False

# 总结. Python 中的集合类型是一种无序容器，不允许有重复运算，由于底层使用了哈希存储，集合中的元素必须是hashable类型。集合与列表最大的区别在于集合中的元素没有顺序、所以不能够通过索引运算访问元素、但是集合可以执行交集、并集、差集等二元运算，也可以通过关系运算符检查两个集合是否存在超集、子集等关系。




"""
Python Learn - Day 13 within 100 Days - Dict
Author: guisongwu
Date  : 2025-10-27
"""
# 迄今为止，我们已经为大家介绍了 Python 中的三种容器型数据类型（列表、元组、集合），但是这些数据类型仍然不足以帮助我们解决所有的问题。例如，我们需要一个变量来保存一个人的多项信息，包括：姓名、年龄、身高、体重、家庭住址、本人手机号、紧急联系人手机号，此时你会发现，我们之前学过的列表、元组和集合类型都不够好使。
# person1 = ['王大锤', 55, 168, 60, '成都市武侯区科华北路62号1栋101', '13122334455', '13800998877'] list 
# person2 = ('王大锤', 55, 168, 60, '成都市武侯区科华北路62号1栋101', '13122334455', '13800998877') tuple
# person3 = {'王大锤', 55, 168, 60, '成都市武侯区科华北路62号1栋101', '13122334455', '13800998877'} set
# 集合肯定是最不合适的，因为集合中不能有重复元素，如果一个人的年龄和体重刚好相同，那么集合中就会少一项信息；同理，如果这个人的手机号和紧急联系人手机号是相同的，那么集合中又会少一项信息。另一方面，虽然列表和元组可以把一个人的所有信息都保存下来，但是当你想要获取这个人的手机号或家庭住址时，你得先知道他的手机号是列表或元组中的第几个元素。总之，在遇到上述的场景时，列表、元组、集合都不是最合适的选择，此时我们需要字典（dictionary）类型，这种数据类型最适合把相关联的信息组装到一起，可以帮助我们解决 Python 程序中为真实事物建模的问题。
# Python 程序中的字典跟现实生活中的字典很像，它以键值对（键和值的组合）的方式把数据组织到一起，我们可以通过键找到与之对应的值并进行操作。就像《新华字典》中，每个字（键）都有与它对应的解释（值）一样，每个字和它的解释合在一起就是字典中的一个条目，而字典中通常包含了很多个这样的条目。
# 创建和使用字典. Python 中创建字典可以使用{}字面量语法，这一点跟上一节课讲的集合是一样的。但是字典的{}中的元素是以键值对的形式存在的，每个元素由:分隔的两个值构成，:前面是键，:后面是值，代码如下所示。
# xinhua = {
#         '麓': '山脚下',
#         '路': '道，往来通行的地方；方面，地区：南～货，外～货；种类：他俩是一～人',
#         '蕗': '甘草的别名',
#         '潞': '潞水，水名，即今山西省的浊漳河；潞江，水名，即云南省的怒江'
#         }
# print(xinhua)
# print(type(xinhua))
# person = {
#         'name': '王大锤',
#         'age': 55,
#         'height': 168,
#         'weight': 60,
#         'addr': '成都市武侯区科华北路62号1栋101', 
#         'tel': '13122334455',
#         'emergence contact': '13800998877'
#         }
# print(person)
# print(type(person))

# 通过上面的代码，相信大家已经看出来了，用字典来保存一个人的信息远远优于使用列表或元组，因为我们可以用:前面的键来表示条目的含义，而:后面就是这个条目所对应的值。当然，如果愿意，我们也可以使用内置函数dict或者是字典的生成式语法来创建字典，代码如下所示。
# dict函数(构造器)中的每一组参数就是字典中的一组键值对
# person = dict(name='王大锤', age=55, height=168, weight=60, addr='成都市武侯区科华北路62号1栋101')
# print(person)  # {'name': '王大锤', 'age': 55, 'height': 168, 'weight': 60, 'addr': '成都市武侯区科华北路62号1栋101'}

# 可以通过Python内置函数zip压缩两个序列并创建字典
# items1 = dict(zip('ABCDE', '12345'))
# print(items1)  # {'A': '1', 'B': '2', 'C': '3', 'D': '4', 'E': '5'}
# items2 = dict(zip('ABCDE', range(1, 10)))
# print(items2)  # {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5}

# 用字典生成式语法创建字典
# items3 = {x: x ** 3 for x in range(1, 6)}
# items4 = {x ** 3 for x in range(1, 6)}
# print(items3)  # {1: 1, 2: 8, 3: 27, 4: 64, 5: 125}
# print(type(items3)) # <class 'dict'>
# print(items4)
# print(type(items4)) # <class 'set'>
# print(items3[4])

# 想知道字典中一共有多少组键值对，仍然是使用len函数；如果想对字典进行遍历，可以用for循环，但是需要注意，for循环只是对字典的键进行了遍历，不过没关系，在学习了字典的索引运算后，我们可以通过字典的键访问它对应的值。
# person = {
#         'name': '王大锤',
#         'age': 55,
#         'height': 168,
#         'weight': 60,
#         'addr': '成都市武侯区科华北路62号1栋101'
#         }
# print(len(person))  # 5
# for key in person:
#     print(key, end=":")
#     print(person[key])

# 字典的运算. 对于字典类型来说，成员运算和索引运算肯定是很重要的，前者可以判定指定的键在不在字典中，后者可以通过键访问对应的值或者向字典中添加新的键值对。值得注意的是，字典的索引不同于列表的索引，列表中的元素因为有属于自己有序号，所以列表的索引是一个整数；字典中因为保存的是键值对，所以字典需要用键去索引对应的值。需要特别提醒大家注意的是，字典中的键必须是不可变类型，例如整数（int）、浮点数（float）、字符串（str）、元组（tuple）等类型，这一点跟集合类型对元素的要求是一样的；很显然，之前我们讲的列表（list）和集合（set）不能作为字典中的键，字典类型本身也不能再作为字典中的键，因为字典也是可变类型，但是列表、集合、字典都可以作为字典中的值，例如：
# person = {
#         'name': '王大锤',
#         'age': 55,
#         'height': 168,
#         'weight': 60,
#         'addr': ['成都市武侯区科华北路62号1栋101', '北京市西城区百万庄大街1号'],
#         'car': {
#             'brand': 'BMW X7',
#             'maxSpeed': '250',
#             'length': 5170,
#             'width': 2000,
#             'height': 1835,
#             'displacement': 3.0
#             }
#         }
# print(person)

# person = {'name': '王大锤', 'age': 55, 'height': 168, 'weight': 60, 'addr': '成都市武侯区科华北路62号1栋101'}
# # 成员运算
# print('name' in person)  # True
# print('tel' in person)   # False
# # 索引运算
# print(person['name'])
# print(person['addr'])
# person['age'] = 25
# person['height'] = 178
# person['tel'] = '13122334455'
# person['signature'] = '你的男朋友是一个盖世垃圾，他会踏着五彩祥云去迎娶你的闺蜜'
# print(person)

# 循环遍历
# person = {'name': '王大锤', 'age': 55, 'height': 168, 'weight': 60, 'addr': '成都市武侯区科华北路62号1栋101'}
# for key in person:
#     print(f'{key}:\t{person[key]}')
# 需要注意，在通过索引运算获取字典中的值时，如指定的键没有在字典中，将会引发KeyError异常。

# 字典的方法. 字典类型的方法基本上都跟字典的键值对操作相关，其中get方法可以通过键来获取对应的值。跟索引运算不同的是，get方法在字典中没有指定的键时不会产生异常，而是返回None或指定的默认值，代码如下所示。
# person = {'name': '王大锤', 'age': 25, 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}
# print(person.get('name'))       # 王大锤
# print(person.get('sex'))        # None
# print(person.get('sex', True))  # True
# 如果需要获取字典中所有的键，可以使用keys方法；如果需要获取字典中所有的值，可以使用values方法。字典还有一个名为items的方法，它会将键和值组装成二元组，通过该方法来遍历字典中的元素也是非常方便的。
# person = {'name': '王大锤', 'age': 25, 'height': 178}
# print(person.keys())    # dict_keys(['name', 'age', 'height'])
# print(type(person.keys())) # <class 'dict_keys'>
# print(person.values())  # dict_values(['王大锤', 25, 178])
# print(person.items())   # dict_items([('name', '王大锤'), ('age', 25), ('height', 178)])
# for key, value in person.items():
#     print(f'{key}:\t{value}')

# 字典的update方法实现两个字典的合并操作。例如，有两个字典x和y，当执行x.update(y)操作时，x跟y相同的键对应的值会被y中的值更新，而y中有但x中没有的键值对会直接添加到x中，代码如下所示。
# person1 = {'name': '王大锤', 'age': 55, 'height': 178}
# person2 = {'age': 25, 'addr': '成都市武侯区科华北路62号1栋101'}
# person1.update(person2)
# print(person1)  # {'name': '王大锤', 'age': 25, 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}
# 如果使用 Python 3.9 及以上的版本，也可以使用|运算符来完成同样的操作，代码如下所示。
# person1 = {'name': '王大锤', 'age': 55, 'height': 178}
# person2 = {'age': 25, 'addr': '成都市武侯区科华北路62号1栋101'}
# person1 |= person2
# print(person1)  # {'name': '王大锤', 'age': 25, 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}

# 可以通过pop或popitem方法从字典中删除元素，前者会返回（获得）键对应的值，但是如果字典中不存在指定的键，会引发KeyError错误；后者在删除元素时，会返回（获得）键和值组成的二元组。字典的clear方法会清空字典中所有的键值对，代码如下所示。
# person = {'name': '王大锤', 'age': 25, 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}
# print(person.pop('age'))  # 25
# print(person)             # {'name': '王大锤', 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}
# print(person.popitem())   # ('addr', '成都市武侯区科华北路62号1栋101')
# print(person)             # {'name': '王大锤', 'height': 178}
# person.clear()
# print(person)             # {}
# 跟列表一样，从字典中删除元素也可以使用del关键字，在删除元素的时候如果指定的键索引不到对应的值，一样会引发KeyError错误，具体的做法如下所示。
# person = {'name': '王大锤', 'age': 25, 'height': 178, 'addr': '成都市武侯区科华北路62号1栋101'}
# del person['age']
# del person['addr']
# print(person)  # {'name': '王大锤', 'height': 178}

# 字典的应用
# 我们通过几个简单的例子来看看如何使用字典类型解决一些实际的问题。
# 例子1：输入一段话，统计每个英文字母出现的次数，按出现次数从高到低输出。
# sentence = input('请输入一段话: ')
# print(type(sentence))
# counter = {} # empty dict
# for ch in sentence:
#     print(ch)
#     if 'A' <= ch <= 'Z' or 'a' <= ch <= 'z':
#         counter[ch] = counter.get(ch, 0) + 1
# sorted_keys = sorted(counter, key=counter.get, reverse=True)
# for key in sorted_keys:
#     print(f'{key} 出现了 {counter[key]} 次.')

# sorted() 是 Python 内置的排序函数，用来对 任意可迭代对象（iterable） 进行排序，返回一个新的已排序列表（不会修改原对象）。
# sorted(iterable, key=None, reverse=False)
# | 参数         | 说明                                  |
# | ---------- | ----------------------------------- |
# | `iterable` | 可迭代对象，比如 list、tuple、str、dict、set 等  |
# | `key`      | 排序时用的函数，用于提取比较关键字                   |
# | `reverse`  | 是否倒序（默认 `False`，即升序；设为 `True` 表示降序） |
# nums = [5, 2, 9, 1]
# print(sorted(nums))
# print(sorted(nums, reverse=True))
# words = ["apple", "pear", "banana", "kiwi"]
# print(sorted(words, key=len))
# scores = {'Alice': 85, 'Bob': 92, 'Charlie': 78}
# print(sorted(scores, key=lambda name: scores[name], reverse=True))
# items = [(1, 3), (2, 1), (4, 2)]
# print(sorted(items, key=lambda x: x[1]))

# 输入：
# Man is distinguished, not only by his reason, but by this singular passion from other animals, which is a lust of the mind, that by a perseverance of delight in the continued and indefatigable generation of knowledge, exceeds the short vehemence of any carnal pleasure.
# 输出：
# e 出现了 27 次.
# n 出现了 21 次.
# a 出现了 18 次.
# i 出现了 18 次.
# s 出现了 16 次.
# t 出现了 16 次.
# o 出现了 14 次.
# h 出现了 13 次.
# r 出现了 10 次.
# d 出现了 9 次.
# l 出现了 9 次.
# g 出现了 6 次.
# u 出现了 6 次.
# f 出现了 6 次.
# c 出现了 6 次.
# y 出现了 5 次.
# b 出现了 5 次.
# m 出现了 4 次.
# p 出现了 3 次.
# w 出现了 2 次.
# v 出现了 2 次.
# M 出现了 1 次.
# k 出现了 1 次.
# x 出现了 1 次.

# 例子2：在一个字典中保存了股票的代码和价格，找出股价大于100元的股票并创建一个新的字典。
# 说明：可以用字典的生成式语法来创建这个新字典。
# stocks = {
#         'AAPL': 191.88,
#         'GOOG': 1186.96,
#         'IBM': 149.24,
#         'ORCL': 48.44,
#         'ACN': 166.89,
#         'FB': 208.09,
#         'SYMC': 21.29
#         }
# stocks2 = {key: value for key, value in stocks.items() if value > 100}
# print(stocks2) # {'AAPL': 191.88, 'GOOG': 1186.96, 'IBM': 149.24, 'ACN': 166.89, 'FB': 208.09}

# 总结. Python 程序中的字典跟现实生活中字典非常像，允许我们以键值对的形式保存数据，再通过键访问对应的值。字典是一种非常有利于数据检索的数据类型，但是需要再次提醒大家，字典中的键必须是不可变类型，列表、集合、字典等类型的数据都不能作为字典的键.




"""
Python Learn - Day 14 within 100 Days - function and modules
Author: guisongwu
Date  : 2025-10-29
"""

# 函数和模块. 在讲解本节课的内容之前，我们先来研究一道数学题，请说出下面的方程有多少组正整数解。
# x1 + x2 + x3 + x4 = 8
# 你可能已经想到了，这个问题其实等同于将 8 个苹果分成四组且每组至少一个苹果有多少种方案，也等价于在分隔 8 个苹果的 7 个间隙之间放入三个隔断将苹果分成四组有多少种方案，所以答案是 C_7^3 = 35 ，其中 C_m^n = \frac{m!}{n!(m-n)!}. 根据之前学习的知识，我们可以用循环做累乘的方式分别计算出 m ! 、 n ! 和 ( m − n  ) ! ，然后再通过除法运算得到组合数 C m n ，代码如下所示。
# m = 7
# n = 3
# # 计算m的阶乘
# fm = 1
# for num in range(1, m + 1):
#     fm *= num
# # 计算n的阶乘
# fn = 1
# for num in range(1, n + 1):
#     fn *= num
# # 计算m-n的阶乘
# fk = 1
# for num in range(1, m - n + 1):
#     fk *= num
# # 计算C(M,N)的值
# print(fm // fn // fk)
# print(fm / fn / fk)

# 不知大家是否注意到，上面的代码中我们做了三次求阶乘的操作，虽然 m 、 n 、 m − n 的值各不相同，但是三段代码并没有实质性的区别，属于重复代码。世界级的编程大师Martin Fowler曾经说过：“代码有很多种坏味道，重复是最坏的一种！”。要写出高质量的代码，首先就要解决重复代码的问题。对于上面的代码来说，我们可以将求阶乘的功能封装到一个称为“函数”的代码块中，在需要计算阶乘的地方，我们只需“调用函数”即可实现对求阶乘功能的复用。
# 定义函数. 数学上的函数通常形如 y = f(x) 或者 z = g(x,y) 这样的形式，在 y=f(x) 中， f 是函数的名字， x 是函数的自变量， y 是函数的因变量；而在 z=g(x,y) 中， g 是函数名， x 和 y 是函数的自变量， z 是函数的因变量。Python 中的函数跟这个结构是一致的，每个函数都有自己的名字、自变量和因变量。我们通常把 Python 函数的自变量称为函数的参数，而因变量称为函数的返回值。
# Python 中可以使用def关键字来定义函数，和变量一样每个函数也应该有一个漂亮的名字，命名规则跟变量的命名规则是一样的（大家赶紧想想我们之前讲过的变量的命名规则）。在函数名后面的圆括号中可以设置函数的参数，也就是我们刚才说的函数的自变量，而函数执行完成后，我们会通过return关键字来返回函数的执行结果，这就是我们刚才说的函数的因变量。如果函数中没有return语句，那么函数会返回代表空值的None。另外，函数也可以没有自变量（参数），但是函数名后面的圆括号是必须有的。一个函数要做的事情（要执行的代码），是通过代码缩进的方式放到函数定义行之后，跟之前分支和循环结构的代码块类似，如下图所示。

# 下面，我们将之前代码中求阶乘的操作放到一个函数中，通过这种方式来重构上面的代码。所谓重构，是在不影响代码执行结果的前提下对代码的结构进行调整，重构之后的代码如下所示。
# 通过关键字def定义求阶乘的函数
# 自变量（参数）num是一个非负整数
# 因变量（返回值）是num的阶乘
# def fac(num):
#     result = 1
#     for n in range(2, num + 1):
#         result *= n
#     return result
# m = 7
# n = 3
# print(fac(m) // fac(n) // fac(m - n))

# 大家可以感受下，上面的代码是不是比之前的版本更加简单优雅。更为重要的是，我们定义的求阶乘函数fac还可以在其他需要求阶乘的代码中重复使用。所以，使用函数可以帮助我们将功能上相对独立且会被重复使用的代码封装起来，当我们需要这些的代码，不是把重复的代码再编写一遍，而是通过调用函数实现对既有代码的复用。事实上，Python 标准库的math模块中，已经有一个名为factorial的函数实现了求阶乘的功能，我们可以直接用import math导入math模块，然后使用math.factorial来调用求阶乘的函数；我们也可以通过from math import factorial直接导入factorial函数来使用它，代码如下所示。
# from math import factorial
# m = 7
# n = 3
# print(factorial(m) // factorial(n) // factorial(m - n))
# 将来我们使用的函数，要么是自定义的函数，要么是 Python 标准库或者三方库中提供的函数，如果已经有现成的可用的函数，我们就没有必要自己去定义，“重复发明轮子”是一件非常糟糕的事情。对于上面的代码，如果你觉得factorial这个名字太长，书写代码的时候不是特别方便，我们在导入函数的时候还可以通过as关键字为其别名。在调用函数的时候，我们可以用函数的别名，而不再使用它之前的名字，代码如下所示。
# from math import factorial as f
# m = 7
# n = 3
# print(f(m) // f(n) // f(m - n))

# 函数的参数
# 位置参数和关键字参数
# 我们再来写一个函数，根据给出的三条边的长度判断是否可以构成三角形，如果可以构成三角形则返回True，否则返回False，代码如下所示。
# 上面make_judgement函数有三个参数，这种参数叫做位置参数，在调用函数时通常按照从左到右的顺序依次传入，而且传入参数的数量必须和定义函数时参数的数量相同，如下所示。
# 如果不想按照从左到右的顺序依次给出a、b、c 三个参数的值，也可以使用关键字参数，通过“参数名=参数值”的形式为函数传入参数，如下所示。
# def make_judgement(a, b, c):
#     """判断三条边的长度能否构成三角形"""
#     return a + b > c and b + c > a and a + c > b
# print(make_judgement(1, 2, 3))  # False
# print(make_judgement(4, 5, 6))  # True
# print(make_judgement(b=2, c=3, a=1))  # False
# print(make_judgement(c=6, b=4, a=5))  # True

# 在定义函数时，我们可以在参数列表中用/设置强制位置参数（positional-only arguments），用*设置命名关键字参数。所谓强制位置参数，就是调用函数时只能按照参数位置来接收参数值的参数；而命名关键字参数只能通过“参数名=参数值”的方式来传递和接收参数，大家可以看看下面的例子。
# /前面的参数是强制位置参数
# def make_judgement(a, b, c, /):
#     """判断三条边的长度能否构成三角形"""
#     return a + b > c and b + c > a and a + c > b
# print(make_judgement(1, 2, 3))
# print(make_judgement(b=2, c=3, a=1))
# 上面的代码会产生TypeError错误，错误信息提示“强制位置参数是不允许给出参数名的”
# TypeError: make_judgement() got some positional-only arguments passed as keyword arguments
# 说明：强制位置参数是 Python 3.8 引入的新特性，在使用低版本的 Python 解释器时需要注意。
# *后面的参数是命名关键字参数
# def make_judgement(*, a, b, c):
#     """判断三条边的长度能否构成三角形"""
#     return a + b > c and b + c > a and a + c > b
# print(make_judgement(b=2, c=3, a=1))
# print(make_judgement(1, 2, 3))
# 上面的代码会产生TypeError错误，错误信息提示“函数没有位置参数但却给了3个位置参数”
# TypeError: make_judgement() takes 0 positional arguments but 3 were given

# 参数的默认值. Python 中允许函数的参数拥有默认值，我们可以把之前讲过的一个例子“CRAPS赌博游戏”（《第07课：分支和循环结构的应用》）中摇色子获得点数的功能封装到函数中，代码如下所示。
# from random import randrange
# 定义摇色子的函数
# 函数的自变量（参数）n表示色子的个数，默认值为2
# 函数的因变量（返回值）表示摇n颗色子得到的点数
# def roll_dice(n=2):
#     total = 0
#     for _ in range(n):
#         total += randrange(1, 7)
#     return total
# print(roll_dice())
# print(roll_dice(3))

# 我们再来看一个更为简单的例子。
# def add(a=0, b=0, c=0):
#     """三个数相加求和"""
#     return a + b + c
# print(add())         # 调用add函数，没有传入参数，那么a、b、c都使用默认值0
# print(add(1))        # 调用add函数，传入一个参数，该参数赋值给变量a, 变量b和c使用默认值0
# print(add(1, 2))     # 调用add函数，传入两个参数，分别赋值给变量a和b，变量c使用默认值0
# print(add(1, 2, 3))  # 调用add函数，传入三个参数，分别赋值给a、b、c三个变量
# 需要注意的是，带默认值的参数必须放在不带默认值的参数之后，否则将产生SyntaxError错误，错误消息是：non-default argument follows default argument，翻译成中文的意思是“没有默认值的参数放在了带默认值的参数后面”。

# 可变参数. Python 语言中可以通过星号表达式语法让函数支持可变参数。所谓可变参数指的是在调用函数时，可以向函数传入0个或任意多个参数。将来我们以团队协作的方式开发商业项目时，很有可能要设计函数给其他人使用，但有的时候我们并不知道函数的调用者会向该函数传入多少个参数，这个时候可变参数就能派上用场。下面的代码演示了如何使用可变位置参数实现对任意多个数求和的add函数，调用函数时传入的参数会保存到一个元组，通过对该元组的遍历，可以获取传入函数的参数。
# 用星号表达式来表示args可以接收0个或任意多个参数
# 调用函数时传入的n个参数会组装成一个n元组赋给args
# 如果一个参数都没有传入，那么args会是一个空元组
# def add(*args):
#     total = 0
#     # 对保存可变参数的元组进行循环遍历
#     for val in args:
#     # 对参数进行了类型检查（数值型的才能求和）
#         if type(val) in (int, float):
#             total += val
#     return total
# print(add())                        # 0
# print(add(1))                       # 1
# print(add(1, 2, 3))                 # 6
# print(add(1, 2, 'hello', 3.45, 6))  # 12.45

# 如果我们希望通过“参数名=参数值”的形式传入若干个参数，具体有多少个参数也是不确定的，我们还可以给函数添加可变关键字参数，把传入的关键字参数组装到一个字典中，代码如下所示。
# 参数列表中的**kwargs可以接收0个或任意多个关键字参数
# 调用函数时传入的关键字参数会组装成一个字典（参数名是字典中的键，参数值是字典中的值）
# 如果一个关键字参数都没有传入，那么kwargs会是一个空字典
# def foo(*args, **kwargs):
#     print(args)
#     print(kwargs)
# foo(3, 2.1, True, name='骆昊', age=43, gpa=4.95)

# 用模块管理函数. 不管用什么样的编程语言来写代码，给变量、函数起名字都是一个让人头疼的问题，因为我们会遇到命名冲突这种尴尬的情况。最简单的场景就是在同一个.py文件中定义了两个同名的函数，如下所示。
# def foo():
#     print('hello, world!')
# def foo():
#     print('goodbye, world!')
# foo()  # goodbye, world!

# 当然上面的这种情况我们很容易就能避免，但是如果项目是团队协作多人开发的时候，团队中可能有多个程序员都定义了名为foo的函数，这种情况下怎么解决命名冲突呢？答案其实很简单，Python 中每个文件就代表了一个模块（module），我们在不同的模块中可以有同名的函数，在使用函数的时候，我们通过import关键字导入指定的模块再使用完全限定名（模块名.函数名）的调用方式，就可以区分到底要使用的是哪个模块中的foo函数，代码如下所示。
# import module1
# import module2
# module1.foo()  # hello, world!
# module2.foo()  # goodbye, world!
# # 在导入模块时，还可以使用as关键字对模块进行别名，这样我们可以使用更为简短的完全限定名。
# import module1 as m1
# import module2 as m2
# m1.foo()  # hello, world!
# m2.foo()  # goodbye, world!
# # 上面两段代码，我们导入的是定义函数的模块，我们也可以使用from...import...语法从模块中直接导入需要使用的函数. 但是，如果我们如果从两个不同的模块中导入了同名的函数，后面导入的函数会替换掉之前的导入，就像下面的代码，调用foo会输出goodbye, world!，因为我们先导入了module1的foo，后导入了module2的foo 。如果两个from...import...反过来写，那就是另外一番光景了。
# from module1 import foo
# from module2 import foo
# foo()  # goodbye, world!
# # 如果想在上面的代码中同时使用来自两个模块的foo函数还是有办法的，大家可能已经猜到了，还是用as关键字对导入的函数进行别名，代码如下所示。
# from module1 import foo as f1
# from module2 import foo as f2
# f1()  # hello, world!
# f2()  # goodbye, world!

# 标准库中的模块和函数. Python 标准库中提供了大量的模块和函数来简化我们的开发工作，我们之前用过的random模块就为我们提供了生成随机数和进行随机抽样的函数；而time模块则提供了和时间操作相关的函数；我们之前用到过的math模块中还包括了计算正弦、余弦、指数、对数等一系列的数学函数。随着我们深入学习 Python 语言，我们还会用到更多的模块和函数。
# Python 标准库中还有一类函数是不需要import就能够直接使用的，我们将其称之为内置函数，这些内置函数不仅有用而且还很常用，下面的表格列出了一部分的内置函数。
# 函数    说明
# abs     返回一个数的绝对值，例如：abs(-1.3)会返回1.3。
# bin     把一个整数转换成以'0b'开头的二进制字符串，例如：bin(123)会返回'0b1111011'。
# chr     将Unicode编码转换成对应的字符，例如：chr(8364)会返回'€'。
# ord     将字符转换成对应的Unicode编码，例如：ord('€')会返回8364。
# hex     将一个整数转换成以'0x'开头的十六进制字符串，例如：hex(123)会返回'0x7b'。
# input   从输入中读取一行，返回读到的字符串。
# len     获取字符串、列表等的长度。
# max     返回多个参数或一个可迭代对象中的最大值，例如：max(12, 95, 37)会返回95。
# min     返回多个参数或一个可迭代对象中的最小值，例如：min(12, 95, 37)会返回12。
# oct     把一个整数转换成以'0o'开头的八进制字符串，例如：oct(123)会返回'0o173'。
# open    打开一个文件并返回文件对象。
# pow     求幂运算，例如：pow(2, 3)会返回8；pow(2, 0.5)会返回1.4142135623730951。
# print   打印输出。
# range   构造一个范围序列，例如：range(100)会产生0到99的整数序列。
# round   按照指定的精度对数值进行四舍五入，例如：round(1.23456, 4)会返回1.2346。
# sum     对一个序列中的项从左到右进行求和运算，例如：sum(range(1, 101))会返回5050。
# type    返回对象的类型，例如：type(10)会返回int；而 type('hello')会返回str。
# 总结. 函数是对功能相对独立且会重复使用的代码的封装。学会使用定义和使用函数，就能够写出更为优质的代码。当然，Python 语言的标准库中已经为我们提供了大量的模块和常用的函数，用好这些模块和函数就能够用更少的代码做更多的事情；如果这些模块和函数不能满足我们的要求，可能就需要自定义函数，然后再通过模块的概念来管理这些自定义函数。



"""
Python Learn - Day 15 within 100 Days - Functions
Author: guisongwu
Date  : 2025-10-30
"""
# 函数应用实战
# 例子1：随机验证码
# 设计一个生成随机验证码的函数，验证码由数字和英文大小写字母构成，长度可以通过参数设置。
# import random
# import string
# ALL_CHARS = string.digits + string.ascii_letters
# a = (random.choice(ALL_CHARS))
# print(a)
# print(type(a)) # <class 'str'>
# b = (random.choices(ALL_CHARS, k=4))
# print(b)
# print(type(b)) # <class 'list'>
# def generate_code(*, code_len=4):
#     """
#     生成指定长度的验证码
#     :param code_len: 验证码的长度(默认4个字符)
#     :return: 由大小写英文字母和数字构成的随机验证码字符串
#     """
#     return ''.join(random.choices(ALL_CHARS, k=code_len))
# 说明1：string模块的digits代表0到9的数字构成的字符串'0123456789'，string模块的ascii_letters代表大小写英文字母构成的字符串'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'。
# 说明2：random模块的sample和choices函数都可以实现随机抽样，sample实现无放回抽样，这意味着抽样取出的元素是不重复的；choices实现有放回抽样，这意味着可能会重复选中某些元素。这两个函数的第一个参数代表抽样的总体，而参数k代表样本容量，需要说明的是choices函数的参数k是一个命名关键字参数，在传参时必须指定参数名。
# for _ in range(5):
#     print(generate_code()) 
# for _ in range(5):
#     print(generate_code(code_len=6))
# 说明：我们设计的generate_code函数的参数是命名关键字参数，由于它有默认值，可以不给它传值，使用默认值4。如果需要给函数传入参数，必须指定参数名code_len。


# 例子2：判断素数. 设计一个判断给定的大于1的正整数是不是质数的函数。质数是只能被1和自身整除的正整数（大于1），如果一个大于 1 的正整数 N 是质数，那就意味着在 2 到 N − 1 之间都没有它的因子。
# def is_prime(num: int) -> bool:
#     """
#     判断一个正整数是不是质数
#     :param num: 大于1的正整数
#     :return: 如果num是质数返回True，否则返回False
#     """
#     for i in range(2, int(num ** 0.5) + 1):
#         if num % i == 0:
#             return False
#     return True
# print(is_prime(10))
# print(is_prime(11))
# print(is_prime(13))
# 说明1：上面is_prime函数的参数num后面的: int用来标注参数的类型，虽然它对代码的执行结果不产生任何影响，但是很好的增强了代码的可读性。同理，参数列表后面的-> bool用来标注函数返回值的类型，它也不会对代码的执行结果产生影响，但是却让我们清楚的知道，调用函数会得到一个布尔值，要么是True，要么是False。
# 说明2：上面的循环并不需要从 2 循环到 N − 1 ，因为如果循环进行到 \sqrt{N} 时，还没有找到 N 的因子，那么 \sqrt{N} 之后也不会出现 N 的因子，大家可以自己想一想这是为什么。


# 例子3：最大公约数和最小公倍数. 设计计算两个正整数最大公约数和最小公倍数的函数。 x 和 y 的最大公约数是能够同时整除 x 和 y 的最大整数，如果 x 和 y 互质，那么它们的最大公约数为 1； x 和 y 的最小公倍数是能够同时被 x 和 y 整除的最小正整数，如果 x 和 y 互质，那么它们的最小公倍数为 x × y 。需要提醒大家注意的是，计算最大公约数和最小公倍数是两个不同的功能，应该设计成两个函数，而不是把两个功能放到同一个函数中。
# def lcm(x: int, y: int) -> int:
#     """求最小公倍数"""
#     return x * y // gcd(x, y)
# def gcd(x: int, y: int) -> int:
#     """求最大公约数"""
#     while y % x != 0:
#         x, y = y % x, x
#     return x
# print(lcm(10, 34))
# 说明：函数之间可以相互调用，上面求最小公倍数的lcm函数调用了求最大公约数的gcd函数，通过 x × y g c d ( x , y  ) 来计算最小公倍数。

# 例子4：数据统计. 假设样本数据保存一个列表中，设计计算样本数据描述性统计信息的函数。描述性统计信息通常包括：算术平均值、中位数、极差（最大值和最小值的差）、方差、标准差、变异系数等，计算公式如下所示。
# 样本均值（sample mean）：
# \bar{x} = \frac{\sum_{i=1}^n x_i}{n}
# 样本方差（sample variance）：
# s^2 = \sum_{i=1}^{n}(x_i-\bar{x})^2 / (n-1)
# 样本标准差（sample standard deviation）：
# s = \sqrt(s^2)
# 变异系数（coefficient of sample variation）：
# CV = s / \bar{x}

def ptp(data):
    """极差（全距）"""
    return max(data) - min(data)
def mean(data):
    """算术平均"""
    return sum(data) / len(data)
def median(data):
    """中位数"""
    temp, size = sorted(data), len(data)
    if size % 2 != 0:
        return temp[size // 2]
    else:
        return mean(temp[size // 2 - 1:size // 2 + 1])
def var(data, ddof=1):
    """方差"""
    x_bar = mean(data)
    temp = [(num - x_bar) ** 2 for num in data]
    return sum(temp) / (len(temp) - ddof)
def std(data, ddof=1):
    """标准差"""
    return var(data, ddof) ** 0.5


def cv(data, ddof=1):
    """变异系数"""
    return std(data, ddof) / mean(data)


def describe(data):
    """输出描述性统计信息"""
    print(f'均值: {mean(data)}')
    print(f'中位数: {median(data)}')
    print(f'极差: {ptp(data)}')
    print(f'方差: {var(data)}')
    print(f'标准差: {std(data)}')
    print(f'变异系数: {cv(data)}')

    说明1：中位数是将数据按照升序或降序排列后位于中间的数，它描述了数据的中等水平。中位数的计算分两种情况：当数据体量$n$为奇数时，中位数是位于 n + 1 2 位置的元素；当数据体量 n 为偶数时，中位数是位于 n 2 和 n 2 + 1 两个位置元素的均值。

    说明2：计算方差和标准差的函数中有一个名为ddof的参数，它代表了可以调整的自由度，默认值为 1。在计算样本方差和样本标准差时，需要进行自由度校正；如果要计算总体方差和总体标准差，可以将ddof参数赋值为 0，即不需要进行自由度校正。

    说明3：describe函数将上面封装好的统计函数组装到一起，用于输出数据的描述性统计信息。事实上，Python 标准库中有一个名为statistics的模块，它已经把获取描述性统计信息的函数封装好了，有兴趣的读者可以自行了解。

例子5：双色球随机选号

我们用函数重构之前讲过的双色球随机选号的例子（《第09课：常用数据结构之列表-2》），将生成随机号码和输出一组号码的功能分别封装到两个函数中，然后通过调用函数实现机选N注号码的功能。

"""
双色球随机选号程序

Author: 骆昊
Version: 1.3
"""
import random

RED_BALLS = [i for i in range(1, 34)]
BLUE_BALLS = [i for i in range(1, 17)]


def choose():
    """
    生成一组随机号码
    :return: 保存随机号码的列表
    """
    selected_balls = random.sample(RED_BALLS, 6)
    selected_balls.sort()
    selected_balls.append(random.choice(BLUE_BALLS))
    return selected_balls


def display(balls):
    """
    格式输出一组号码
    :param balls: 保存随机号码的列表
    """
    for ball in balls[:-1]:
        print(f'\033[031m{ball:0>2d}\033[0m', end=' ')
    print(f'\033[034m{balls[-1]:0>2d}\033[0m')


n = int(input('生成几注号码: '))
for _ in range(n):
    display(choose())

    说明：大家看看display(choose())这行代码，这里我们先通过choose函数获得一组随机号码，然后把choose函数的返回值作为display函数的参数，通过display函数将选中的随机号码显示出来。重构之后的代码逻辑非常清晰，代码的可读性更强了。如果有人为你封装了这两个函数，你仅仅是函数的调用者，其实你根本不用关心choose函数和display函数的内部实现，你只需要知道调用choose函数可以生成一组随机号码，而调用display函数传入一个列表，就可以输出这组号码。将来我们使用各种各样的 Python 三方库时，我们也根本不关注它们的底层实现，我们需要知道的仅仅是调用哪个函数可以解决问题。

总结
在写代码尤其是开发商业项目的时候，一定要有意识的将相对独立且重复使用的功能封装成函数，这样不管是自己还是团队的其他成员都可以通过调用函数的方式来使用这些功能，减少工作中那些重复且乏味的劳动。










