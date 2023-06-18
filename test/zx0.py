from tkinter import N


def staff(wo='原神芭芭拉',
          url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0',
          nu=90):
    from zx_h import (change, downLoadjpg, gethtmlbypyppteer,
                      gethtmlbyrequests, gethtmlbyselenium)

    # se = input_set(se)
    url = change(wo, url)
    html = gethtmlbypyppteer(url)
    # html= gethtmlbyrequests(url)
    # html = gethtmlbyselenium(url)
    print(html)
    downLoadjpg(wo, nu, html)


def img2(pa='img/game.png', spa='img/game2.png'):
    from PIL import Image
    im = Image.open(pa)
    # im = im.convert('RGB')
    im = im.convert('RGBA')
    im.save(spa, quality=100)


def panch():
    from time import sleep

    from pyautogui import doubleClick, hotkey, keyDown, keyUp, typewrite
    url = r'http://my.lzu.edu.cn:8080/login?service=http://my.lzu.edu.cn'
    path = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    keyDown('win')
    hotkey('r')
    keyUp('win')
    typewrite(url)
    hotkey('enter')
    sleep(0.5)
    doubleClick(1346, 440)
    sleep(0.5)
    doubleClick(1346, 510)
    sleep(0.5)
    doubleClick(1346, 510)
    sleep(0.5)
    doubleClick(282, 573)
    sleep(1)
    doubleClick(1675, 687)
    sleep(0.5)
    doubleClick(1672, 1000)
    sleep(0.5)
    doubleClick(1920, 20)


def news():
    import re

    import requests
    print('国际新闻')
    url = 'https://news.china.com/international/'
    r = requests.get(url)
    r.encoding = 'utf-8'
    html = r.text
    alist = re.findall(
        r'.html" target="_blank">.+</a></h3><span class="item_info">', html)
    blist = re.findall(r'<em class="item_source">.+</em>', html)
    for i in range(20):
        print('({})'.format(i + 1), end='')
        print(blist[i][24:-42], end='于')
        print(blist[i][-15:-5], end='报道：')
        print(alist[i][23:-33])
    print('兰大新闻')
    url = u'http://jwc.lzu.edu.cn/'
    r = requests.get(url)
    r.encoding = 'GB2312'
    html = r.text
    li = re.findall(
        r'<li><a title=\'(.+?)\'\s\shref=\".*?/lzupage(.+?)\"', html)
    for i in range(-10, 0):
        print('({})'.format(i + 11), end='')
        print(li[i][0] + '...http://jwc.lzu.edu.cn/lzupage' + li[i][1])
    input('input for leave')


def getpotion():
    import time

    import pyautogui

    try:
        while True:
            x, y = pyautogui.position()
            rgb = pyautogui.screenshot().getpixel((x, y))
            posi = 'x:' + str(x).rjust(4) + ' y:' + \
                   str(y).rjust(4) + '  RGB:' + str(rgb)
            print('\r', posi, end='')
            time.sleep(0.5)

    except KeyboardInterrupt:
        print('已退出！')


def timer():
    from os import startfile
    from time import localtime, sleep
    while True:
        if localtime().tm_hour == 7 and localtime().tm_min == 25:
            startfile(r"D:\Programming\User\dist\panch.exe")
            sleep(50)
            startfile(
                r"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\腾讯软件\QQ音乐\QQ音乐.lnk")
            sleep(5)
            startfile(r"D:\Programming\User\dist\news.exe")
            sleep(5)
            startfile(
                r"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\OneNote.lnk")
            sleep(60)
        sleep(40)


def fileStatisticsbyPython_zx():
    from os import listdir, path

    from pandas import DataFrame

    pa = str(input('postion(\\):'))
    li = []
    for i in listdir(pa):
        if path.isdir(pa + '\\' + i):
            li.append(i)
    lc = []
    for i in li:
        lt = list('=HYPERLINK(\"' + pa + '\\' + i + '\\' + j +
                  '\",\"' + j + '\")' for j in listdir(pa + '\\' + i))
        lt.insert(0, str(len(listdir(pa + '\\' + i))))
        lc.append(lt)
    df = DataFrame(data=lc, index=li)
    df.to_excel(pa + '\\' + 'fileStatisticsbyPython_zx.xlsx', na_rep=' ')
    print(df)
    input('input for leave')


def install():
    import os
    libs = ["pyppeteer", "pdf2docx", "numpy", "matplotlib", "pillow", "sklearn", "requests", "jieba", "pyperclip",
            "wheel", "pyautogui", "sympy", "pyinstaller", "Cpython", "selenium", "paddleocr", "opencv-python -i",
            "pandas",
            "wordcloud", "ipython", "paddlepaddle", "openpyxl"]
    # libs=['pyautogui']

    try:
        for lib in libs:
            os.system("python -m pip install " + lib)
            print("Successful")
    except BaseException:
        print("Failed Somehow")


def temtimer():
    from os import startfile
    from time import sleep
    startfile(r"D:\Programming\User\dist\panch.exe")
    sleep(50)
    startfile(
        r"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\腾讯软件\QQ音乐\QQ音乐.lnk")
    sleep(5)
    startfile(r"D:\Programming\User\dist\news.exe")
    sleep(5)
    startfile(r"C:\ProgramData\Microsoft\Windows\Start Menu\Programs\OneNote.lnk")


def originalMaze(n=20, m=30):
    import numpy as np
    m0 = np.random.rand(n, m)
    for i in range(n):
        for j in range(m):
            if m0[i][j] < 0.25:
                m0[i][j] = 0
            else:
                m0[i][j] = 1
    m0[n // 2][m // 2] = 1
    return m0


def mazeSolution(n=20, m=30):
    import random

    import numpy as np
    m1 = np.ones((n, m))
    p = [n // 2, m // 2]
    while True:
        try:
            t = random.randint(1, 4)
            if p[0] == 1:
                return m1
            elif p[1] == 1:
                return m1
            elif t % 4 == 0:
                p[0] -= 2
                m1[p[0] + 1][p[1]] = 0
                m1[p[0]][p[1]] = 0
            elif t % 4 == 1:
                p[0] += 2
                m1[p[0] - 1][p[1]] = 0
                m1[p[0]][p[1]] = 0
            elif t % 4 == 2:
                p[1] -= 2
                m1[p[0]][p[1] + 1] = 0
                m1[p[0]][p[1]] = 0
            elif t % 4 == 3:
                p[1] += 2
                m1[p[0]][p[1] - 1] = 0
                m1[p[0]][p[1]] = 0
        except IndexError:
            m1[n // 2][m // 2] = 2
            return m1


def prin(mn, n, m):
    for i in range(n):
        print(str(list(mn[i])).replace(".", '.').replace(
            '0', ' ').replace('1', '#').replace(',', '').replace('2', '@'))


def mazemain(m=40, n=40):
    m0 = originalMaze(m=m, n=n)
    m1 = mazeSolution(m=m, n=n)
    m2 = m0 * m1
    prin(mn=m2, m=m, n=n)
    print(
        "-------------------------------------------------------------------------------------------------------------------------------------------")
    prin(mn=m1, m=m, n=n)


def isprime(n=100):
    import numpy as np
    is_prime = np.ones((n,), dtype=bool)
    is_prime[:2] = 0
    N_max = int(np.sqrt(len(is_prime)))
    for j in range(2, N_max):
        is_prime[2 * j::j] = False
    return is_prime

    # for i in range(len(ni)):
    #     f += ni[-i]*x**(i)
    # plot((a, b), f)

    # return f


def lineStas(x, y):
    import matplotlib.pyplot as plt
    from scipy import stats
    from seaborn import set
    set()
    plt.scatter(x, y)
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
    getmodel = list(map(slope * x + intercept, x))
    plt.plot(x, getmodel)


def polyFit(x=[], y=[], deg=1, text='', xlabel='', ylabel=''):
    import matplotlib.pyplot as plt
    import pandas as pd
    from numpy import corrcoef, poly1d, polyfit
    from seaborn import set
    set()
    plt.scatter(x, y)
    model = poly1d(polyfit(x, y, deg))
    r = pd.DataFrame(corrcoef(x, y))
    plt.plot(x, model(x))
    plt.xlabel(xlabel=xlabel)
    plt.ylabel(ylabel=ylabel)
    plt.title(str(model)+"\n{}\n".format(text))
    # plt.title(str(model) +
    #           '\nCorrelation coefficient \n(take rows in order as new rows, first-order correlation, multi-order ignore)\n' + str(
    #     r))
    plt.show()
    return model, r


# import inspect
#
# for name, obj in inspect.getmembers(inspect):
#     if inspect.isclass(obj):
#         print(name, str(obj))
#     if inspect.isfunction(obj):
# print(name, str(obj))
if __name__ == '__main__':
    # a = [40, 45, 50, 55, 60]
    # b = [0.39, 0.44, 0.49, 0.54, 0.58]
    # polyFit(a, b, 1)
    # import numpy as np
    # news()
    # news()
    # panch()
    # import matplotlib.pyplot as plt
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    # # 支持中文
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    # from matplotlib import rc
    # font = {'family': 'SimHei', "size": 24}
    # rc('font', **font)  # 一次定义终身使用
    # from gethtml.gethtmlbypyppeteer import gethtmlbypyppteer
    # f = hehe(n=1)
    # x = np.linspace(0.5, 3.5, 7, endpoint=True)*1e-3*9.793
    # y = np.array([12.5, 25.1, 38.4, 51.5, 66.6, 79.5, 91.5])
    # deg = 1
    # print(x, y)
    # Polyfit(x, y, deg)
    # staff(wo='利兹与青鸟', nu=30, url='https://yandex.com/images/search?family=yes&rpt=imageview&cbir_page=similar&url=https%3A%2F%2Fyandex-images.clstorage.net%2FU5l21Ro55%2F5e2d48nVvkh%2Ftv4qHX8gnMkCD0xATFyjmg29CMMkKQSyaWePtHYmuaXOl-F7d5rZoEm4y4V-BD3tLx4KY6AUcIJLytsJKJA5rrinbY21b9FAVh2dR1zeY8vKMoSKbSeY-iSo1b_85eqkDNbmqbOKWi8ELFucKNygOLYAi6wrjj8JWwLayGER8S94rjNK_GiagkBJCxkCXObETvbABi3q1nzib1ywfSUh7QiXLDp6iNsjB3LVnKtd6zczhQQlrYT5wN0BG80l0XCn9-H6RXFuEQUBiQqbBBRjxYq1y8ZuZBF7pCDZsL_p62pMW-SxsklYLQ890thlUq7_soyEJenEconSyBWCqZd9ZDdh-sc-4NlIBgqD0sbbLcEHcE_C62Rc-Opq0fOxYmjrQ5mr8_BEGihNfR3UpZ3leL5AiyzlCbWD0AWUQmfX_ig6Y_2Cs-WRAUCBgR-JUOWHA7TIgCyinfVkYxfzMqhv40HfLLY1RRBpjXUUGqXRIz3xBU2tLQR_iFnFlMVsmvpuMq0wgz6p0ARPgspaTN7uyQ_3A4Sg4l18qurYcHWjLmxBnCBzv4uY7oA02JogV24_Mo5NI6CLuIqZwNwK41L8ZHLqdM37oJWBB8qGlkZY40OJMUZHLqSdOqXmH3b6oi-izN2iP7KBWKgLMNzQZ9PqObFJROJhTbMLGQ8ciGGY-ez7bHTCdSFdRkwKANZKmK5JSjTJQqvqlfzsr5h2sS_h5kDepjmyylUmz7-dkCHa5PW9D8ptK8S_yBlPX4CkmP0mN6D3xfKvUY2ODwoeQtityEV5RwBvpp3yLKYXPfWh5ysHHiE3sk1fLc73EJBhWWN9947PKuFL-oKUzVACIVq4oXkiMgA66FQMSs3L0sTYJoxJ_A5GZC6SdCFu1nj6K26uSVUm-zrFFCoCeJJVplig_3BHyCqiTXCKmcGcCqRbvmG_7bNF_eccRgjFhNuBl2lAQH2Pwuvlmzqop1c5dc&cbir_id=1495804%2FgbGpJxKlDatVo2JaWcaXdg1998&crop=0.15%3B0.15%3B0.85%3B0.85')
    # news()
    # n=1000000
    # print(isprime(1000))
    # li=isprime(n)
    # for i in range(n):
    #     if li[i] :
    #         print(i)
    # print(gethtmlbypyppteer())
    # news()
    # mazemain(n=80, m=80
    # news()
    # from matplotlib import pyplot as plt
    # import numpy as np

    # x1 = np.linspace(1, 10, 20)
    # print(x1)
    # y1 = x1 * x1 + 2
    # fig = plt.figure()
    # axes = fig.add_axes([0.1, 0.1, 0.9, 0.9])
    # axes.plot(x1, y1, 'r')
    # plt.show()

    # def fileStatisticsbyPython_zx(pah='./'):
    #     from os import listdir
    #     from os.path import isdir, join, split
    #     from pandas import DataFrame
    #     lit = []
    #     for i in listdir(pah):
    #         if isdir(join(pah, i)):
    #             lit.append(join(pah, i))
    #     print(lit)
    #     lc = []
    #     for i in lit:
    #         lt = list('=HYPERLINK(\"' + join(pah, i, j) +
    #                   '\",\"' + j + '\")' for j in listdir(join(pah, i)))
    #         lt.insert(0, str(len(listdir(join(pah, i)))))
    #         lc.append(lt)
    #     df = DataFrame(data=lc, index=lit)
    #     df.to_excel(join(pah, 'fileStatisticsbyPython_zx.xlsx'), na_rep=' ')
    #     print(df)

    # def findAllfiles(pah):
    #     from os import walk, path
    #     lit = []
    #     for foldName, subfolders, filenames in walk(pah):
    #         for filename in filenames:
    #             try:
    #                 lit.append(path.join(foldName, filename))
    #             except Exception as e:
    #                 print(e)
    #                 pass
    #     return lit

    # def fileSort(pah=[], want_lit=['pdf', 'pptx'], file_lit=[]):
    #     from os import mkdir, rename
    #     from os.path import join, split
    #     for i in want_lit:
    #         try:
    #             i = mkdir(join(pah, i))
    #         except Exception as e:
    #             print(e)
    #         finally:
    #             for j in file_lit:
    #                 (_, j_) = split(j)
    #                 if i in j_:
    #                     rename(j, join(pah, i, j_))

    # if __name__ == '__main__':
    #     pah = '/home/zhangxin/桌面'
    #     # absolute path
    #     lit = findAllfiles(pah=pah)
    #     fileSort(file_lit=lit, pah=pah)
    #     fileStatisticsbyPython_zx(pah=pah)

    # # 稳定双共轭梯度法
    # import numpy as np
    # n = 4
    # converged = np.array([False])
    # v = np.zeros((n, 1))
    # r = np.zeros((n, 1))
    # r0_hat = np.zeros((n, 1))
    # p = np.zeros((n, 1))
    # s = np.zeros((n, 1))
    # t = np.zeros((n, 1))
    # A = np.random.rand(n, n)
    # b = np.random.rand(n, 1)
    # x = np.ones((n, 1))
    # tol = 1e-6
    # limit = 1e6
    # print('系数矩阵 A:\n', A)
    # print('右边向量 b:\n', b)
    # print('初始猜测值 x0:\n', x)
    # r[:] = b[:]-np.dot(A, x)
    # r0_hat[:] = r[:]
    # rho0 = alpha = w = v[:] = p[:] = 1.0
    # rho1 = np.dot(np.transpose(r0_hat), r)
    # iters = 0
    # while(True):
    #     iters = iters+1
    #     converged[:] = (np.linalg.norm(r) < tol*np.linalg.norm(b))
    #     if converged == True or iters == limit:
    #         break
    #     # if iters < 5:
    #     #     print("先验值x({}):\n{}".format(iters, x))
    #     beta = (rho1/rho0)*(alpha/w)
    #     p[:, 0] = r[:, 0]+beta*(p[:, 0]-w*v[:, 0])
    #     v[:] = np.dot(A, p)
    #     alpha = rho1/np.dot(np.transpose(r0_hat), v)
    #     s[:, 0] = r[:, 0]-alpha*v[:, 0]
    #     t[:] = np.dot(A, s)
    #     w = np.dot(np.transpose(t), s)/np.dot(np.transpose(t), t)
    #     rho0 = rho1
    #     rho1 = -w*np.dot(np.transpose(r0_hat), t)
    #     x[:, 0] = x[:, 0]+alpha*p[:, 0]+w*s[:, 0]
    #     r[:, 0] = s[:, 0]-w*t[:, 0]

    # print('迭代次数 iters:\n', iters)
    # print('解向量 x:\n', x)
    # import numpy as np
    # n = 4
    # converged = np.array([False])
    # v = np.zeros((n, 1))
    # r = np.zeros((n, 1))
    # r0_hat = np.zeros((n, 1))
    # p = np.zeros((n, 1))
    # s = np.zeros((n, 1))
    # t = np.zeros((n, 1))
    # A = np.random.rand(n, n)
    # b = np.random.rand(n, 1)
    # x = np.ones((n, 1))
    # tol = 1e-6
    # limit = 1e6
    # print('系数矩阵 A:\n', A)
    # print('右边向量 b:\n', b)
    # print('初始猜测值 x0:\n', x)
    # r[:] = b[:]-np.dot(A, x)
    # r0_hat[:] = r[:]
    # rho0 = alpha = w = v[:] = p[:] = 1.0
    # rho1 = np.dot(np.transpose(r0_hat), r)
    # iters = 0
    # while(True):
    #     iters = iters+1
    #     converged[:] = (np.linalg.norm(r) < tol*np.linalg.norm(b))
    #     if converged == True or iters == limit:
    #         break
    #     # if iters < 5:
    #     #     print("先验值x({}):\n{}".format(iters, x))
    #     beta = (rho1/rho0)*(alpha/w)
    #     p[:, 0] = r[:, 0]+beta*(p[:, 0]-w*v[:, 0])
    #     v[:] = np.dot(A, p)
    #     alpha = rho1/np.dot(np.transpose(r0_hat), v)
    #     s[:, 0] = r[:, 0]-alpha*v[:, 0]
    #     t[:] = np.dot(A, s)
    #     w = np.dot(np.transpose(t), s)/np.dot(np.transpose(t), t)
    #     rho0 = rho1
    #     rho1 = -w*np.dot(np.transpose(r0_hat), t)
    #     x[:, 0] = x[:, 0]+alpha*p[:, 0]+w*s[:, 0]
    #     r[:, 0] = s[:, 0]-w*t[:, 0]

    # print('迭代次数 iters:\n', iters)
    # print('解向量 x:\n', x)

    # # -*- coding: utf-8 -*-

    # import os
    # from pdfminer.pdfparser import PDFParser, PDFDocument
    # from pdfminer.pdfinterp import PDFTextExtractionNotAllowed
    # from pdfminer.layout import *
    # from pdfminer.converter import PDFPageAggregator
    # from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
    # import time
    # import random
    # import hashlib
    # import sys
    # import importlib
    # importlib.reload(sys)

    # # **********翻译部分********************
    # def fanyi(query):
    #     import http.client
    #     import hashlib
    #     import urllib
    #     import random
    #     import json
    #     appid = ''  # !!!!补充
    #     secretKey = ''  # !!!!补充
    #     httpClient = None
    #     myurl = '/api/trans/vip/translate'
    #     q = query
    #     fromLang = 'auto'
    #     toLang = 'zh'
    #     salt = random.randint(32768, 65536)

    #     sign = appid + q + str(salt) + secretKey
    #     m1 = hashlib.md5()
    #     m1.update(sign.encode())
    #     sign = m1.hexdigest()
    #     myurl = myurl + '?appid=' + appid + '&q=' + urllib.parse.quote(
    #         q) + '&from=' + fromLang + '&to=' + toLang + '&salt=' + str(salt) + '&sign=' + sign
    #     # print(urllib.parse.quote(q))
    #     try:
    #         httpClient = http.client.HTTPConnection('api.fanyi.baidu.com')
    #         httpClient.request('GET', myurl)
    #         # response是HTTPResponse对象
    #         response = httpClient.getresponse()
    #         html = response.read()  # bytes
    #         # print("html:  ",type(html),html)
    #         html_str = html.decode()  # bytes to str
    #         # print("html_str:  ",type(html_str),html_str)
    #         html_dict = json.loads(html_str)  # str to dict
    #         # print("html_dist:  ",type(html_dict),html_str)
    #         # result_ori = html_dict["trans_result"][0]["src"]
    #         # result_tar = html_dict["trans_result"][0]["dst"]
    #         # print(html_dict["trans_result"])
    #         result_tar = ''
    #         for i in html_dict["trans_result"]:
    #             result_tar += i["dst"]
    #         # print(result_ori, " --> ", result_tar)
    #         print("翻译文本: " + result_tar)
    #         print("*" * 100)
    #         return result_tar
    #     except Exception as e:
    #         print(e)
    #         return ''
    #     finally:
    #         if httpClient:
    #             httpClient.close()

    # '''
    # 解析pdf文件，获取文件中包含的各种对象
    # '''

    # # 解析pdf文件函数

    # def parse(pdf_path):
    #     textName = pdf_path.split('\\')[-1].split('.')[0] + '.txt'
    #     fp = open(pdf_path, 'rb')  # 以二进制读模式打开
    #     # 用文件对象来创建一个pdf文档分析器
    #     parser = PDFParser(fp)
    #     # 创建一个PDF文档
    #     doc = PDFDocument()
    #     # 连接分析器 与文档对象
    #     parser.set_document(doc)
    #     doc.set_parser(parser)

    #     # 提供初始化密码
    #     # 如果没有密码 就创建一个空的字符串
    #     doc.initialize()

    #     # 检测文档是否提供txt转换，不提供就忽略
    #     if not doc.is_extractable:
    #         raise PDFTextExtractionNotAllowed
    #     else:
    #         # 创建PDf 资源管理器 来管理共享资源
    #         rsrcmgr = PDFResourceManager()
    #         # 创建一个PDF设备对象
    #         laparams = LAParams()
    #         device = PDFPageAggregator(rsrcmgr, laparams=laparams)
    #         # 创建一个PDF解释器对象
    #         interpreter = PDFPageInterpreter(rsrcmgr, device)

    #         # 用来计数页面，图片，曲线，figure，水平文本框等对象的数量
    #         num_page, num_image, num_curve, num_figure, num_TextBoxHorizontal = 0, 0, 0, 0, 0

    #         # 循环遍历列表，每次处理一个page的内容
    #         for page in doc.get_pages():  # doc.get_pages() 获取page列表
    #             num_page += 1  # 页面增一
    #             print("\r\n>> 当前页：", num_page)
    #             interpreter.process_page(page)
    #             # 接受该页面的LTPage对象
    #             layout = device.get_result()
    #             for x in layout:
    #                 if isinstance(x, LTImage):  # 图片对象
    #                     num_image += 1
    #                 if isinstance(x, LTCurve):  # 曲线对象
    #                     num_curve += 1
    #                 if isinstance(x, LTFigure):  # figure对象
    #                     num_figure += 1
    #                 if isinstance(x, LTTextBoxHorizontal):  # 获取文本内容
    #                     num_TextBoxHorizontal += 1  # 水平文本框对象增一
    #                     results = x.get_text()
    #                     print(results.replace('\n', ''))
    #                     # 保存文本内容
    #                     with open(textName, 'a+', encoding='utf8') as f:
    #                         results = x.get_text()
    #                         f.write(results.replace('\n', '') + '\n')
    #         print('对象数量：\n', '页面数：%s\n' % num_page, '图片数：%s\n' % num_image, '曲线数：%s\n' % num_curve, '水平文本框：%s\n'
    #               % num_TextBoxHorizontal)

    # if __name__ == '__main__':

    #     pdf_path = r'sympy-docs-pdf-1.10.1.pdf'
    #     rootPath = '\\'.join(pdf_path.split('\\')[:-1]) if "\\" in pdf_path else ''
    #     textName = pdf_path.split('\\')[-1].split('.')[0] + '.txt'
    #     print(">> 当前文件：", os.path.join(rootPath, textName))
    #     if os.path.exists(os.path.join(rootPath, textName)):
    #         print(">> 删除：", textName)
    #         os.remove(os.path.join(rootPath, textName))
    #     if os.path.exists(os.path.join(rootPath, "translate.txt")):
    #         print(">> 删除：", "translate.txt")
    #         os.remove(os.path.join(rootPath, "translate.txt"))

    #     parse(pdf_path)

    #     with open(textName, 'r', encoding='utf8') as f:
    #         content = f.read()
    #         results = content.split('.')
    #         for i in results:
    #             res = fanyi(i)
    #             with open("translate.txt", 'a+', encoding='utf8') as fp:
    #                 fp.write(res + '\n')

    #             time.sleep(1)

    # from iog_reader import iog_read
    # from os import listdir
    # from re import match
    # # import jax.numpy as np
    # import numpy as np
    # import pandas as pd
    # import gvar as gv
    # import lsqfit
    # import matplotlib.pyplot as plt

    # pd.set_option("display.max_rows", None)

    # # @ti.kernel
    # def get_c2pt(
    #     iog_path="Data",
    #     intrptr=["cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"],
    #     N=24,
    # ):
    #     iog_files = [
    #         iog_path + "/" + file
    #         for file in listdir(iog_path)
    #         if match(r".*dat.iog", file) != None
    #     ]
    #     df_list = []
    #     for iog_file in iog_files:
    #         df0 = iog_read(iog_file, intrptr)
    #         # print(
    #         #     df0[
    #         #         (df0["tsrc"] == "505046")
    #         #         & (df0["lnk"] == "9200505")
    #         #         & (df0["snksmr"] == "1")
    #         #     ][["Re", "cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"]]
    #         # )
    #         df1 = df0[
    #             (df0["tsrc"] == "505046")
    #             & (df0["lnk"] == "9200505")
    #             & (df0["snksmr"] == "1")
    #         ]["Re"]
    #         df1 = df1.reset_index(drop=True)
    #         df2 = np.array([df1[i::N] for i in range(N)])
    #         df3 = np.mean(df2, axis=1)
    #         df_list.append(df3)
    #     print(np.array(df_list).shape)
    #     return np.array(df_list)

    # def jcknf_c2pt(c2pt=None):
    #     Ncnfg = c2pt.shape[0]
    #     c2pt_jcknf = (np.sum(c2pt, axis=0) - c2pt) / (Ncnfg - 1)
    #     return c2pt_jcknf

    # def make_mdls(t_dctnry, p):
    #     mdls = {}
    #     ts = t_dctnry["c2pt"]
    #     mdls["c2pt"] = np.exp(p["n0"] + p["n1"] * ts)
    #     return mdls

    # def fit_c2tp(c2pt_jcknf, make_mdls):
    #     Ncnfg = c2pt_jcknf.shape[0]
    #     T = c2pt_jcknf.shape[1]
    #     T_start = 3
    #     T_end = 10
    #     t_ary = np.array(range(T_start, T_end))
    #     c2pt_jcknf_avg = c2pt_jcknf[:, T_start:T_end]
    #     c2pt_avg_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
    #     c2pt_avg_cov = (Ncnfg - 1) * np.cov(np.transpose(c2pt_jcknf_avg, axes=(1, 0)))
    #     tsep_dctnry = {"c2pt": t_ary}
    #     c2pt_dctnry = {"c2pt": gv.gvar(c2pt_avg_cntrl, c2pt_avg_cov)}

    #     p0 = {
    #         "n0": -25,
    #         "n1": -3,
    #     }
    #     fit = lsqfit.nonlinear_fit(data=(tsep_dctnry, c2pt_dctnry), fcn=make_mdls, p0=p0)
    #     print(fit.format(True))
    # #     t_ary = fit.data[0]["c2pt"]
    # #     c2pt_avg_dat_gvar = fit.data[1]["c2pt"]
    # #     c2pt_avg_dat_cntrl = np.array([c2.mean for c2 in c2pt_avg_dat_gvar])
    # #     c2pt_avg_dat_err = np.array([c2.sdev for c2 in c2pt_avg_dat_gvar])
    # #     t_lst = np.linspace(10, T - 10, 50)
    # #     c2pt_fit_fcn_gvar = fit.fcn({"c2pt": t_lst}, fit.p)["c2pt"]
    # #     c2pt_fit_fcn_cntrl = np.array([c2.mean for c2 in c2pt_fit_fcn_gvar])
    # #     c2pt_fit_fcn_err = np.array([c2.sdev for c2 in c2pt_fit_fcn_gvar])
    # #     c2pt_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
    # #     c2pt_err = np.sqrt(Ncnfg - 1) * np.std(c2pt_jcknf_avg, axis=0)

    # #     plt.figure(dpi=600)
    # #     plt.errorbar(
    # #         t_ary,
    # #         c2pt_cntrl,
    # #         yerr=c2pt_err,
    # #         fmt="bo",
    # #         label="$C_2$",
    # #     )
    # #     plt.errorbar(
    # #         t_ary,
    # #         c2pt_avg_dat_cntrl,
    # #         yerr=c2pt_avg_dat_err,
    # #         fmt="go",
    # #         label="frwrd/bckwrd avg. $C_2$",
    # #     )
    # #     plt.plot(t_lst, c2pt_fit_fcn_cntrl, color="b", label="best fit")
    # #     plt.fill_between(
    # #         t_lst,
    # #         c2pt_fit_fcn_cntrl - c2pt_fit_fcn_err,
    # #         c2pt_fit_fcn_cntrl + c2pt_fit_fcn_err,
    # #     )
    # #     plt.xlabel("t/a")
    # #     plt.ylabel("$C_2$")
    # #     plt.legend(
    # #         loc="upper center",
    # #         frameon=True,
    # #         fancybox=True,
    # #         markerfirst=True,
    # #     )
    # #     plt.savefig("c2pt_fit.png")
    # #     plt.show()
    # def main():
    #     # c2pt=get_c2pt(iog_path='/public/home/sunpeng/chen_c/22_04_13_run_chroma/L24x72_m0235_3pt_mom4_st_20_sm_link/data')
    #     c2pt = get_c2pt(iog_path="Data")
    #     c2pt_jcknf = jcknf_c2pt(c2pt=c2pt)
    #     fit_c2tp(c2pt_jcknf=c2pt_jcknf, make_mdls=make_mdls)

    # # %time throws =main()
    # if __name__ == "__main__":
    #     main()
    # # -*- coding: utf-8 -*-
    # """
    # Created on Fri May  6 14:51:56 2022

    # @author: 15149
    # """
    # # import numpy as np

    # from numpy import poly1d, polyfit, corrcoef, array

    # def Polyfit(x, y, deg=1):
    #     import matplotlib.pyplot as plt
    #     # from seaborn import set
    #     # import pandas as pd

    #     # set()
    #     # 以上为预设，非必要
    #     plt.scatter(x, y)
    #     model = poly1d(polyfit(x, y, deg))
    #     # r = pd.DataFrame(corrcoef(x, y))
    #     r = corrcoef(x, y)
    #     plt.plot(x, model(x))
    #     plt.title(str(model) +
    #               '\nCorrelation coefficient \n(take rows in order as new rows, first-order correlation, multi-order ignore)\n'+str(r))
    #     plt.show()
    #     return model, r

    # if __name__ == '__main__':

    #     # import matplotlib.pyplot as plt

    #     # # 支持中文
    #     # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    #     # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    #     # from matplotlib import rc
    #     # font = {'family': 'SimHei', "size": 24}
    #     # rc('font', **font)  # 一次定义终身使用
    #     # from gethtml.gethtmlbypyppeteer import gethtmlbypyppteer
    #     # f = hehe(n=1)
    #     print('输入样式：\n[0.0048965,0.009793,0.0146895,0.019586,0.0244825,0.029379,0.0342755]\n[12.5,25.1,38.4,51.5,66.6,79.5,91.5]\n1\n')
    #     x = input('x\n')
    #     y = input('y\n')
    #     deg = input('deg(拟合阶数)\n')
    #     x = array(eval(x))
    #     y = array(eval(y))
    #     deg = int(deg)
    #     # print(x, y, deg)
    #     print(Polyfit(x, y, deg))
    #     input("input for out")
    # # -*- coding: utf-8 -*-
    # """
    # Created on Fri May  6 14:51:56 2022

    # @author: 15149
    # """
    # import numpy as np

    # def Polyfit(x, y, deg=1):
    #     import matplotlib.pyplot as plt
    #     from seaborn import set
    #     import pandas as pd
    #     from numpy import poly1d, polyfit, corrcoef
    #     set()
    #     # 以上为预设，非必要
    #     plt.scatter(x, y)
    #     model = poly1d(polyfit(x, y, deg))
    #     r = pd.DataFrame(corrcoef(x, y))
    #     plt.plot(x, model(x))
    #     plt.title(str(model) +
    #               '\nCorrelation coefficient \n(take rows in order as new rows, first-order correlation, multi-order ignore)\n'+str(r))
    #     plt.show()
    #     return model, r

    # if __name__ == '__main__':

    #     # import matplotlib.pyplot as plt

    #     # # 支持中文
    #     # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    #     # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    #     # from matplotlib import rc
    #     # font = {'family': 'SimHei', "size": 24}
    #     # rc('font', **font)  # 一次定义终身使用
    #     # from gethtml.gethtmlbypyppeteer import gethtmlbypyppteer
    #     # f = hehe(n=1)
    #     # print('输入样式：\n[0.0048965,0.009793,0.0146895,0.019586,0.0244825,0.029379,0.0342755]\n[12.5,25.1,38.4,51.5,66.6,79.5,91.5]\n1\n')
    #     # x = input('x\n')
    #     # y = input('y\n')
    #     # deg = input('deg(拟合阶数)\n')
    #     y = np.array([40.17, 26.78, 20.09, 16.09, 13.39, 11.48])/100.0

    #     x = np.array([1.146, 0.509, 0.333, 0.206, 0.137, 0.0979])
    #     print('x:', x, '\ny:', y)
    #     x = np.log(x)
    #     y = np.log(y)
    #     print('ln(x)', x, '\nln(y)', y)
    #     # deg = int(deg)
    #     deg = 1
    #     # print(x, y, deg, type(x), type(y), type(deg))
    #     print(Polyfit(x, y, deg))
    #     # input("input for out")
    # [0.0048965, 0.009793, 0.0146895, 0.019586, 0.0244825, 0.029379, 0.0342755]

    # def hehe(x1=[0.9460666, 0.11884469, -0.48295271, -0.93720798, -1.19210023],
    #          y1=[-8.87466906, -1.28013417, -1.56781624, -1.79095979, -1.97328135], y2=[-0.87466906, -1.28013417, -1.56781624, -1.79095979, -1.97328135], n=2):
    #     from numpy import polyfit, array
    #     from sympy import symbols, lambdify, ln
    #     import seaborn as sns
    #     import matplotlib.pyplot as plt
    #     sns.set()
    #     plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    #     plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    #     x = symbols('x')
    #     f1 = symbols('f1', functions=True)
    #     f2 = symbols('f2', functions=True)
    #     f = lambdify(x, ln(x))
    #     # a = [0.2, 0.092, 0.051, 0.032, 0.022, 0.016]
    #     # b = linspace(76.6/2, 76.6/7, 6, endpoint=True)/100

    #     # print(a, b)
    #     # a, b = f(array(a)*10), f(array(b))
    #     x1, y1, y2 = array(x1), array(y1), array(y2)
    #     n1 = polyfit(x1, y1, n)
    #     n2 = polyfit(x1, y2, n)

    #     # print(a, b, ni)
    #     f1 = 0
    #     for i in range(len(n1)):
    #         f1 += n1[-(i + 1)] * x ** (i)
    #     F1 = lambdify(x, f1)
    #     f2 = 0
    #     for i in range(len(n2)):
    #         f2 += n2[-(i + 1)] * x ** (i)
    #     F2 = lambdify(x, f2)
    #     plt.plot(x1, y1, 'r:o', label='升温实验')
    #     plt.plot(x1, F1(x1), 'r-o', label='升温拟合')
    #     plt.plot(x1, y2, 'b:o', label='降温实验')
    #     plt.plot(x1, F2(x1), 'b-o', label='降温拟合')
    #     plt.title('升温拟合：f1={}\n降温拟合：f2={}'.format(f1, f2))
    #     plt.xlabel('温度值（C。）')
    #     plt.ylabel('长度值（mm-3）')
    #     plt.legend()
    #     plt.show()
    #     return F1, F2

    # import math

    # def inputSet(a):
    #     from pyautogui import prompt, alert
    #     for ai in a:
    #         ay = a[ai]
    #         t = type(a[ai])
    #         while True:
    #             a[ai] = ay
    #             try:
    #                 p = prompt(
    #                     "{}|默认:{};类型:{};修改:输入'@'开始修改;退出:输入'#'准备退出".format(ai, a[ai], t))
    #                 if p == '':
    #                     break
    #                 if p == '#':
    #                     exit(0)
    #                 if p == '@':
    #                     a[ai] = prompt(ai)
    #                     if not isinstance(eval(a[ai]), t):
    #                         alert('输入有误')
    #                     else:
    #                         alert('输入完成')
    #                         break
    #             except Exception as e:
    #                 alert(e)
    #                 alert('重新输入')
    #     return a

    # def downLoadwebp(wo='原神芭芭拉', nu=30, html=''):
    #     from re import findall
    #     from os import mkdir
    #     from io import BytesIO
    #     from PIL import Image
    #     from requests import get
    #     from getpass import getuser
    #     from pyautogui import alert
    #     try:
    #         mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    #     except FileExistsError:
    #         pass
    #     li = findall(r'\"(https://[^;]+?)\"', html)
    #     li = list(set(li))
    #     alert('已找到{}张图片'.format(len(li)))
    #     for i in range(nu):
    #         try:
    #             pb = BytesIO(get(li[i]).content)
    #             print(li[i])
    #             im = Image.open(pb)
    #             if im.mode == "RGBA":
    #                 im.load()  # required for png.split()
    #                 background = Image.new("RGB", im.size, (255, 255, 255))
    #                 background.paste(im, mask=im.split()[3])
    #             path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(
    #                 getuser(), wo, wo, i + 1, nu)
    #             im.save(path, 'JPEG')
    #         except Exception as e:
    #             print(e)
    #             continue
    #     alert('下载完毕(可能漏下)')

    # def getHtmlbyRequests(
    #         url=r'https://image.baidu.com/search/index?tn=baiduimage&ipn=r&ct=201326592&cl=2&lm=-1&st=-1&sf=1&fmq=&pv=&ic=0&nc=1&z=&se=1&showtab=0&fb=0&width=&height=&face=0&istype=2&ie=utf-8&fm=index&pos=history&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89',
    #         endcoding='utf-8'):
    #     from requests import get
    #     headers = {
    #         'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
    #         'Accept-Encoding': 'gzip, deflate, br',
    #         'Accept-Language': 'zh-CN,zh;q=0.9,en-US;q=0.8,en;q=0.7',
    #         'Cache-Control': 'max-age=0',
    #         'Connection': 'keep-alive',
    #         'sec-ch-ua': '"Google Chrome";v="89", "Chromium";v="89", ";Not A Brand";v="99"',
    #         'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/89.0.4389.82 Safari/537.36'}
    #     params = {
    #         'kd': '张三'
    #     }
    #     rep = get(url, headers=headers, params=params)
    #     rep.encoding = endcoding
    #     return rep.text

    # def findAllfiles(pa="D:\\Monitor"):
    #     from os import walk, path
    #     f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'w')
    #     for foldName, subfolders, filenames in walk(pa):
    #         for filename in filenames:
    #             try:
    #                 f.write(path.join(foldName, filename) + '\r')
    #             except Exception as e:
    #                 print(e)
    #     f.close()
    #     print('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt')

    # def findallfiles(pa="D:\\Monitor"):
    #     from os import walk, path
    #     li = []
    #     for foldName, subfolders, filenames in walk(pa):
    #         for filename in filenames:
    #             try:
    #                 li.append(path.join(foldName, filename))
    #             except Exception as e:
    #                 print(e)
    #                 pass
    #     return li

    # def pdftoword(fileName):
    #     import re
    #     from pdf2docx import Converter
    #     pdf_file = fileName
    #     name = re.findall(r'(.*?)\.', pdf_file)[0]
    #     docx_file = f'{name}.docx'

    #     cv = Converter(pdf_file)
    #     cv.convert(docx_file, start=0, end=None)
    #     cv.close()

    # def getPotion():
    #     import pyautogui
    #     import time

    #     try:
    #         while True:
    #             x, y = pyautogui.position()
    #             rgb = pyautogui.screenshot().getpixel((x, y))
    #             posi = 'x:' + str(x).rjust(4) + ' y:' + \
    #                    str(y).rjust(4) + '  RGB:' + str(rgb)
    #             print('\r', posi, end='')
    #             time.sleep(0.5)

    #     except KeyboardInterrupt:
    #         print('已退出！')

    #     # downloadjpg(se['wo'], se['nu'], html)

    # def gethtmlbyselenium(
    #         url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    #     from selenium import webdriver
    #     from time import sleep
    #     # options = Options()
    #     # options.add_argument('--headless')
    #     # options.add_argument('--disable-gpu')
    #     # driver = webdriver.Edge(options=options)
    #     driver = webdriver.Edge()
    #     driver.get(url)
    #     for i in range(5):
    #         js = "window.scrollTo(0,5000);"
    #         driver.execute_script(js)
    #         sleep(0.5)
    #     html = driver.page_source
    #     driver.close()
    #     driver.quit()
    #     return html

    # def downLoadjpg(wo='原神芭芭拉', nu=30, html=''):
    #     from re import findall
    #     from os import mkdir
    #     from requests import get
    #     from getpass import getuser
    #     from pyautogui import alert
    #     try:
    #         mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    #     except FileExistsError:
    #         pass
    #     li = findall(r'https://[^;]+?\.[pngje]{3,4}', html)
    #     alert('已找到{}张图片'.format(len(li)))
    #     for i in range(nu):
    #         try:
    #             pb = get(li[i]).content
    #             path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(
    #                 getuser(), wo, wo, i + 1, nu)
    #             f = open(path, 'wb')
    #             f.write(pb)
    #             f.close()
    #         except Exception as e:
    #             print(e)
    #             continue
    #     alert('下载完毕(可能漏下)')

    # def gethtmlbypyppteer(
    #         url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    #     import asyncio
    #     from pyppeteer import launch
    #     from time import sleep

    #     async def asgethtml(url):
    #         browser = await launch({
    #             'headless': False,
    #             'devtools': True,
    #             'args': [
    #                 '--disable-extensions',
    #                 '--hide-scrollbars',
    #                 '--disable-bundled-ppapi-flash',
    #                 '--mute-audio',
    #                 '--no-sandbox',
    #                 '--disable-setuid-sandbox',
    #                 '--disable-gpu',
    #             ],
    #             'dumpio': True,
    #         })
    #         page = await browser.newPage()
    #         await page.setUserAgent(
    #             "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.99 Safari/537.36 Edg/97.0.1072.69"
    #         )
    #         await page.evaluateOnNewDocument('()=>{Object.defineProperties(navigator,{webdriver:{get:()=>false}})}')
    #         await page.goto(url)
    #         for i in range(5):
    #             await page.evaluate('window.scrollTo(0,5000)')
    #             sleep(0.5)
    #         text = await page.content()
    #         await browser.close()
    #         return text

    #     m = asyncio.ensure_future(asgethtml(url))
    #     asyncio.get_event_loop().run_until_complete(m)
    #     return m.result()

    # def fileStatisticsbyPython_zx():
    #     from os import listdir, path

    #     from pandas import DataFrame

    #     # pa = str(input('postion(\\):'))
    #     pa = r'D:\Monitor'
    #     # pa = r"D:\Monitor - 副本"
    #     li = []
    #     for i in listdir(pa):
    #         if path.isdir(pa + '\\' + i):
    #             li.append(i)
    #     lc = []
    #     for i in li:
    #         lt = list('=HYPERLINK(\"' + pa + '\\' + i + '\\' + j + '\",\"' + j + '\")' for j in
    #                   listdir(pa + '\\' + i))
    #         for k in lt:
    #             print(k)
    #             if k[-5:-2] == 'cfg':
    #                 print('ok')
    #                 lt.remove(k)
    #         # lt.insert(0, str(len(listdir(pa + '\\' + i))))
    #         lt.insert(0, str(len(lt)))
    #         lc.append(lt)
    #     df = DataFrame(data=lc, index=li)
    #     df.to_excel(pa + '\\' + 'fileStatisticsbyPython_zx.xlsx', na_rep=' ')
    #     print(df)
    #     input('input for leave')

    # def change(wo='原神',
    #            url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0'):
    #     return url.replace('%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89', wo)

    # def ha(pa=r"D:\Programming\User\txt\city.txt"):
    #     import hashlib

    #     m = hashlib.md5()
    #     f = open(pa, 'rb')
    #     m.update(f.read())
    #     f.close()
    #     md5_value = m.hexdigest()
    #     # print(md5_value)
    #     return md5_value

    # def fileTree(pa=r"D:\Programming\User", de=0):
    #     from os.path import isdir
    #     from os import listdir
    #     if de == 0:
    #         print('文件夹:' + pa)
    #     for i in listdir(pa):
    #         print('|    ' * de + '+--' + i)
    #         di = pa + '/' + i
    #         if isdir(di):
    #             fileTree(di, de + 1)

    # def filetree(pa=r"D:\User", de=0, spa=r"D:\Programming\User"):
    #     def fi(pa, de):
    #         from os.path import isdir
    #         from os import listdir
    #         if de == 0:
    #             f.write('文件夹:' + pa + '\n')
    #         for i in listdir(pa):
    #             f.write('|    ' * de + '+--' + i + '\n')
    #             di = pa + '/' + i
    #             if isdir(di):
    #                 fi(di, de + 1)

    #     f = open(spa + '\\' + '_'.join(pa.replace(':',
    #              '').split(sep='\\')) + '_filetree2.txt', 'w')
    #     fi(pa, de)
    #     f.close()

    # def Gravitation(m1=1, m2=1, r=1, G=6.67e-11):
    #     return G * m1 * m2 / (r ** 2)

    # def ElectricForce(e1=1, e2=1, r=1, k=8.99e9):
    #     return k * e1 * e2 / (r ** 2)

    # def SurroundSpeed(m=1, r=1, f=1):
    #     return (f * r / m) ** 0.5

    # def ElectricFieldStrength(q=1, r=1, k=8.99e9):
    #     return q * k / (r * r)

    # def ElectronSpeed(u=1, e=1.6e-19, m=9.1e-31):
    #     return (2 * u * e / m) ** 0.5

    # if __name__ == '__main__':
    #     # print(gethtmlbypyppteer())
    #     print(math.log(2.7))
    #     print('{:e}'.format(ElectronSpeed(u=300)))

    # if __name__ == '__main__':
    #     hehe()
