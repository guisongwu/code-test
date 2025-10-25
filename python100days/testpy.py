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




























