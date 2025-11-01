# ==========================
# Test for Python Programmer
# ==========================

# ---------------------------- color -----------------------------
# a = 100
# print(f"\033[030m a={a}\033[0m") # 30 for black
# print(f"\033[031m a={a}\033[0m") # 31 for red
# print(f"\033[032m a={a}\033[0m") # 32 for green
# print(f"\033[033m a={a}\033[0m") # 33 for yellow
# print(f"\033[034m a={a}\033[0m") # 34 for blue
# print(f"\033[035m a={a}\033[0m") # 35 for purple
# print(f"\033[036m a={a}\033[0m") # 36 for cyan
# print(f"\033[037m a={a}\033[0m") # 37 for white

# --------------------------- randrange --------------------------
# import random
# # a random integer in [0, 10)
# print(random.randrange(10))
# print(random.randrange(0, 10))
# print(random.randrange(0, 10, 1))
# # a random integer in [0, 10) with step 2, maybe 0 2 4 6 8
# print(random.randrange(0, 10, 2))

# ---------------------------- join ------------------------------
# ----- str.join(iterable)
# words = ['Hello', 'World', 'Python']
# result = ' '.join(words)
# print(result)  # Hello World Python
# result = '-'.join(words)
# print(result)  # Hello-World-Python
# result = ', '.join(words)
# print(result)  # Hello, World, Python
# numbers = [1, 2, 3, 4] # 连接数字列表（需要先转换为字符串）
# result1 = '-'.join(map(str, numbers))
# result2 = '-'.join(str(i) for i in numbers)
# print(result1)  # 1-2-3-4
# print(result2)  # 1-2-3-4
# fruits = ('apple', 'banana', 'cherry')
# result = ' and '.join(fruits)
# print(result)  # apple and banana and cherry

# -------------------------- type int ------------------------------
# import sys
# def show_int_sizes():
#     numbers = [0, 1, 100, 1000, 10**10, 10**100]
#     for num in numbers:
#         print(f"{num}: {sys.getsizeof(num)} bytes")
# show_int_sizes()






