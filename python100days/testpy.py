# ==========================
# Test for Python Programmer
# ==========================

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






