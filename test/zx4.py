from iog_reader import iog_read
from os import listdir
from re import match
# import jax.numpy as np
import numpy as np
import pandas as pd
import gvar as gv
import lsqfit
import matplotlib.pyplot as plt

pd.set_option("display.max_rows", None)

# @ti.kernel
def get_c2pt(
    iog_path="Data",
    intrptr=["cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"],
    N=24,
):
    iog_files = [
        iog_path + "/" + file
        for file in listdir(iog_path)
        if match(r".*dat.iog", file) != None
    ]
    df_list = []
    for iog_file in iog_files:
        df0 = iog_read(iog_file, intrptr)
        # print(
        #     df0[
        #         (df0["tsrc"] == "505046")
        #         & (df0["lnk"] == "9200505")
        #         & (df0["snksmr"] == "1")
        #     ][["Re", "cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"]]
        # )
        df1 = df0[
            (df0["tsrc"] == "505046")
            & (df0["lnk"] == "9200505")
            & (df0["snksmr"] == "1")
        ]["Re"]
        df1 = df1.reset_index(drop=True)
        df2 = np.array([df1[i::N] for i in range(N)])
        df3 = np.mean(df2, axis=1)
        df_list.append(df3)
    print(np.array(df_list).shape)
    return np.array(df_list)


def jcknf_c2pt(c2pt=None):
    Ncnfg = c2pt.shape[0]
    c2pt_jcknf = (np.sum(c2pt, axis=0) - c2pt) / (Ncnfg - 1)
    return c2pt_jcknf


def make_mdls(t_dctnry, p):
    mdls = {}
    ts = t_dctnry["c2pt"]
    mdls["c2pt"] = np.exp(p["n0"] + p["n1"] * ts)
    return mdls


def fit_c2tp(c2pt_jcknf, make_mdls):
    Ncnfg = c2pt_jcknf.shape[0]
    T = c2pt_jcknf.shape[1]
    T_start = 3
    T_end = 10
    t_ary = np.array(range(T_start, T_end))
    c2pt_jcknf_avg = c2pt_jcknf[:, T_start:T_end]
    c2pt_avg_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
    c2pt_avg_cov = (Ncnfg - 1) * np.cov(np.transpose(c2pt_jcknf_avg, axes=(1, 0)))
    tsep_dctnry = {"c2pt": t_ary}
    c2pt_dctnry = {"c2pt": gv.gvar(c2pt_avg_cntrl, c2pt_avg_cov)}

    p0 = {
        "n0": -25,
        "n1": -3,
    }
    fit = lsqfit.nonlinear_fit(data=(tsep_dctnry, c2pt_dctnry), fcn=make_mdls, p0=p0)
    print(fit.format(True))
#     t_ary = fit.data[0]["c2pt"]
#     c2pt_avg_dat_gvar = fit.data[1]["c2pt"]
#     c2pt_avg_dat_cntrl = np.array([c2.mean for c2 in c2pt_avg_dat_gvar])
#     c2pt_avg_dat_err = np.array([c2.sdev for c2 in c2pt_avg_dat_gvar])
#     t_lst = np.linspace(10, T - 10, 50)
#     c2pt_fit_fcn_gvar = fit.fcn({"c2pt": t_lst}, fit.p)["c2pt"]
#     c2pt_fit_fcn_cntrl = np.array([c2.mean for c2 in c2pt_fit_fcn_gvar])
#     c2pt_fit_fcn_err = np.array([c2.sdev for c2 in c2pt_fit_fcn_gvar])
#     c2pt_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
#     c2pt_err = np.sqrt(Ncnfg - 1) * np.std(c2pt_jcknf_avg, axis=0)


#     plt.figure(dpi=600)
#     plt.errorbar(
#         t_ary,
#         c2pt_cntrl,
#         yerr=c2pt_err,
#         fmt="bo",
#         label="$C_2$",
#     )
#     plt.errorbar(
#         t_ary,
#         c2pt_avg_dat_cntrl,
#         yerr=c2pt_avg_dat_err,
#         fmt="go",
#         label="frwrd/bckwrd avg. $C_2$",
#     )
#     plt.plot(t_lst, c2pt_fit_fcn_cntrl, color="b", label="best fit")
#     plt.fill_between(
#         t_lst,
#         c2pt_fit_fcn_cntrl - c2pt_fit_fcn_err,
#         c2pt_fit_fcn_cntrl + c2pt_fit_fcn_err,
#     )
#     plt.xlabel("t/a")
#     plt.ylabel("$C_2$")
#     plt.legend(
#         loc="upper center",
#         frameon=True,
#         fancybox=True,
#         markerfirst=True,
#     )
#     plt.savefig("c2pt_fit.png")
#     plt.show()
def main():
    # c2pt=get_c2pt(iog_path='/public/home/sunpeng/chen_c/22_04_13_run_chroma/L24x72_m0235_3pt_mom4_st_20_sm_link/data')
    c2pt = get_c2pt(iog_path="Data")
    c2pt_jcknf = jcknf_c2pt(c2pt=c2pt)
    fit_c2tp(c2pt_jcknf=c2pt_jcknf, make_mdls=make_mdls)


# %time throws =main()
if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:51:56 2022

@author: 15149
"""
# import numpy as np

from numpy import poly1d, polyfit, corrcoef, array


def Polyfit(x, y, deg=1):
    import matplotlib.pyplot as plt
    # from seaborn import set
    # import pandas as pd

    # set()
    # 以上为预设，非必要
    plt.scatter(x, y)
    model = poly1d(polyfit(x, y, deg))
    # r = pd.DataFrame(corrcoef(x, y))
    r = corrcoef(x, y)
    plt.plot(x, model(x))
    plt.title(str(model) +
              '\nCorrelation coefficient \n(take rows in order as new rows, first-order correlation, multi-order ignore)\n'+str(r))
    plt.show()
    return model, r


if __name__ == '__main__':

    # import matplotlib.pyplot as plt

    # # 支持中文
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    # from matplotlib import rc
    # font = {'family': 'SimHei', "size": 24}
    # rc('font', **font)  # 一次定义终身使用
    # from gethtml.gethtmlbypyppeteer import gethtmlbypyppteer
    # f = hehe(n=1)
    print('输入样式：\n[0.0048965,0.009793,0.0146895,0.019586,0.0244825,0.029379,0.0342755]\n[12.5,25.1,38.4,51.5,66.6,79.5,91.5]\n1\n')
    x = input('x\n')
    y = input('y\n')
    deg = input('deg(拟合阶数)\n')
    x = array(eval(x))
    y = array(eval(y))
    deg = int(deg)
    # print(x, y, deg)
    print(Polyfit(x, y, deg))
    input("input for out")
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 14:51:56 2022

@author: 15149
"""
import numpy as np


def Polyfit(x, y, deg=1):
    import matplotlib.pyplot as plt
    from seaborn import set
    import pandas as pd
    from numpy import poly1d, polyfit, corrcoef
    set()
    # 以上为预设，非必要
    plt.scatter(x, y)
    model = poly1d(polyfit(x, y, deg))
    r = pd.DataFrame(corrcoef(x, y))
    plt.plot(x, model(x))
    plt.title(str(model) +
              '\nCorrelation coefficient \n(take rows in order as new rows, first-order correlation, multi-order ignore)\n'+str(r))
    plt.show()
    return model, r


if __name__ == '__main__':

    # import matplotlib.pyplot as plt

    # # 支持中文
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    # from matplotlib import rc
    # font = {'family': 'SimHei', "size": 24}
    # rc('font', **font)  # 一次定义终身使用
    # from gethtml.gethtmlbypyppeteer import gethtmlbypyppteer
    # f = hehe(n=1)
    # print('输入样式：\n[0.0048965,0.009793,0.0146895,0.019586,0.0244825,0.029379,0.0342755]\n[12.5,25.1,38.4,51.5,66.6,79.5,91.5]\n1\n')
    # x = input('x\n')
    # y = input('y\n')
    # deg = input('deg(拟合阶数)\n')
    y = np.array([40.17, 26.78, 20.09, 16.09, 13.39, 11.48])/100.0

    x = np.array([1.146, 0.509, 0.333, 0.206, 0.137, 0.0979])
    print('x:', x, '\ny:', y)
    x = np.log(x)
    y = np.log(y)
    print('ln(x)', x, '\nln(y)', y)
    # deg = int(deg)
    deg = 1
    # print(x, y, deg, type(x), type(y), type(deg))
    print(Polyfit(x, y, deg))
    # input("input for out")
[0.0048965, 0.009793, 0.0146895, 0.019586, 0.0244825, 0.029379, 0.0342755]

def hehe(x1=[0.9460666, 0.11884469, -0.48295271, -0.93720798, -1.19210023],
         y1=[-8.87466906, -1.28013417, -1.56781624, -1.79095979, -1.97328135], y2=[-0.87466906, -1.28013417, -1.56781624, -1.79095979, -1.97328135], n=2):
    from numpy import polyfit, array
    from sympy import symbols, lambdify, ln
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set()
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    x = symbols('x')
    f1 = symbols('f1', functions=True)
    f2 = symbols('f2', functions=True)
    f = lambdify(x, ln(x))
    # a = [0.2, 0.092, 0.051, 0.032, 0.022, 0.016]
    # b = linspace(76.6/2, 76.6/7, 6, endpoint=True)/100

    # print(a, b)
    # a, b = f(array(a)*10), f(array(b))
    x1, y1, y2 = array(x1), array(y1), array(y2)
    n1 = polyfit(x1, y1, n)
    n2 = polyfit(x1, y2, n)

    # print(a, b, ni)
    f1 = 0
    for i in range(len(n1)):
        f1 += n1[-(i + 1)] * x ** (i)
    F1 = lambdify(x, f1)
    f2 = 0
    for i in range(len(n2)):
        f2 += n2[-(i + 1)] * x ** (i)
    F2 = lambdify(x, f2)
    plt.plot(x1, y1, 'r:o', label='升温实验')
    plt.plot(x1, F1(x1), 'r-o', label='升温拟合')
    plt.plot(x1, y2, 'b:o', label='降温实验')
    plt.plot(x1, F2(x1), 'b-o', label='降温拟合')
    plt.title('升温拟合：f1={}\n降温拟合：f2={}'.format(f1, f2))
    plt.xlabel('温度值（C。）')
    plt.ylabel('长度值（mm-3）')
    plt.legend()
    plt.show()
    return F1, F2

import math


def inputSet(a):
    from pyautogui import prompt, alert
    for ai in a:
        ay = a[ai]
        t = type(a[ai])
        while True:
            a[ai] = ay
            try:
                p = prompt(
                    "{}|默认:{};类型:{};修改:输入'@'开始修改;退出:输入'#'准备退出".format(ai, a[ai], t))
                if p == '':
                    break
                if p == '#':
                    exit(0)
                if p == '@':
                    a[ai] = prompt(ai)
                    if not isinstance(eval(a[ai]), t):
                        alert('输入有误')
                    else:
                        alert('输入完成')
                        break
            except Exception as e:
                alert(e)
                alert('重新输入')
    return a


def downLoadwebp(wo='原神芭芭拉', nu=30, html=''):
    from re import findall
    from os import mkdir
    from io import BytesIO
    from PIL import Image
    from requests import get
    from getpass import getuser
    from pyautogui import alert
    try:
        mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    except FileExistsError:
        pass
    li = findall(r'\"(https://[^;]+?)\"', html)
    li = list(set(li))
    alert('已找到{}张图片'.format(len(li)))
    for i in range(nu):
        try:
            pb = BytesIO(get(li[i]).content)
            print(li[i])
            im = Image.open(pb)
            if im.mode == "RGBA":
                im.load()  # required for png.split()
                background = Image.new("RGB", im.size, (255, 255, 255))
                background.paste(im, mask=im.split()[3])
            path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(
                getuser(), wo, wo, i + 1, nu)
            im.save(path, 'JPEG')
        except Exception as e:
            print(e)
            continue
    alert('下载完毕(可能漏下)')


def getHtmlbyRequests(
        url=r'https://image.baidu.com/search/index?tn=baiduimage&ipn=r&ct=201326592&cl=2&lm=-1&st=-1&sf=1&fmq=&pv=&ic=0&nc=1&z=&se=1&showtab=0&fb=0&width=&height=&face=0&istype=2&ie=utf-8&fm=index&pos=history&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89',
        endcoding='utf-8'):
    from requests import get
    headers = {
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'zh-CN,zh;q=0.9,en-US;q=0.8,en;q=0.7',
        'Cache-Control': 'max-age=0',
        'Connection': 'keep-alive',
        'sec-ch-ua': '"Google Chrome";v="89", "Chromium";v="89", ";Not A Brand";v="99"',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/89.0.4389.82 Safari/537.36'}
    params = {
        'kd': '张三'
    }
    rep = get(url, headers=headers, params=params)
    rep.encoding = endcoding
    return rep.text


def findAllfiles(pa="D:\\Monitor"):
    from os import walk, path
    f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'w')
    for foldName, subfolders, filenames in walk(pa):
        for filename in filenames:
            try:
                f.write(path.join(foldName, filename) + '\r')
            except Exception as e:
                print(e)
    f.close()
    print('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt')


def findallfiles(pa="D:\\Monitor"):
    from os import walk, path
    li = []
    for foldName, subfolders, filenames in walk(pa):
        for filename in filenames:
            try:
                li.append(path.join(foldName, filename))
            except Exception as e:
                print(e)
                pass
    return li


def pdftoword(fileName):
    import re
    from pdf2docx import Converter
    pdf_file = fileName
    name = re.findall(r'(.*?)\.', pdf_file)[0]
    docx_file = f'{name}.docx'

    cv = Converter(pdf_file)
    cv.convert(docx_file, start=0, end=None)
    cv.close()


def getPotion():
    import pyautogui
    import time

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

    # downloadjpg(se['wo'], se['nu'], html)


def gethtmlbyselenium(
        url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    from selenium import webdriver
    from time import sleep
    # options = Options()
    # options.add_argument('--headless')
    # options.add_argument('--disable-gpu')
    # driver = webdriver.Edge(options=options)
    driver = webdriver.Edge()
    driver.get(url)
    for i in range(5):
        js = "window.scrollTo(0,5000);"
        driver.execute_script(js)
        sleep(0.5)
    html = driver.page_source
    driver.close()
    driver.quit()
    return html


def downLoadjpg(wo='原神芭芭拉', nu=30, html=''):
    from re import findall
    from os import mkdir
    from requests import get
    from getpass import getuser
    from pyautogui import alert
    try:
        mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    except FileExistsError:
        pass
    li = findall(r'https://[^;]+?\.[pngje]{3,4}', html)
    alert('已找到{}张图片'.format(len(li)))
    for i in range(nu):
        try:
            pb = get(li[i]).content
            path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(
                getuser(), wo, wo, i + 1, nu)
            f = open(path, 'wb')
            f.write(pb)
            f.close()
        except Exception as e:
            print(e)
            continue
    alert('下载完毕(可能漏下)')


def gethtmlbypyppteer(
        url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    import asyncio
    from pyppeteer import launch
    from time import sleep

    async def asgethtml(url):
        browser = await launch({
            'headless': False,
            'devtools': True,
            'args': [
                '--disable-extensions',
                '--hide-scrollbars',
                '--disable-bundled-ppapi-flash',
                '--mute-audio',
                '--no-sandbox',
                '--disable-setuid-sandbox',
                '--disable-gpu',
            ],
            'dumpio': True,
        })
        page = await browser.newPage()
        await page.setUserAgent(
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.99 Safari/537.36 Edg/97.0.1072.69"
        )
        await page.evaluateOnNewDocument('()=>{Object.defineProperties(navigator,{webdriver:{get:()=>false}})}')
        await page.goto(url)
        for i in range(5):
            await page.evaluate('window.scrollTo(0,5000)')
            sleep(0.5)
        text = await page.content()
        await browser.close()
        return text

    m = asyncio.ensure_future(asgethtml(url))
    asyncio.get_event_loop().run_until_complete(m)
    return m.result()


def fileStatisticsbyPython_zx():
    from os import listdir, path

    from pandas import DataFrame

    # pa = str(input('postion(\\):'))
    pa = r'D:\Monitor'
    # pa = r"D:\Monitor - 副本"
    li = []
    for i in listdir(pa):
        if path.isdir(pa + '\\' + i):
            li.append(i)
    lc = []
    for i in li:
        lt = list('=HYPERLINK(\"' + pa + '\\' + i + '\\' + j + '\",\"' + j + '\")' for j in
                  listdir(pa + '\\' + i))
        for k in lt:
            print(k)
            if k[-5:-2] == 'cfg':
                print('ok')
                lt.remove(k)
        # lt.insert(0, str(len(listdir(pa + '\\' + i))))
        lt.insert(0, str(len(lt)))
        lc.append(lt)
    df = DataFrame(data=lc, index=li)
    df.to_excel(pa + '\\' + 'fileStatisticsbyPython_zx.xlsx', na_rep=' ')
    print(df)
    input('input for leave')


def change(wo='原神',
           url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0'):
    return url.replace('%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89', wo)


def ha(pa=r"D:\Programming\User\txt\city.txt"):
    import hashlib

    m = hashlib.md5()
    f = open(pa, 'rb')
    m.update(f.read())
    f.close()
    md5_value = m.hexdigest()
    # print(md5_value)
    return md5_value


def fileTree(pa=r"D:\Programming\User", de=0):
    from os.path import isdir
    from os import listdir
    if de == 0:
        print('文件夹:' + pa)
    for i in listdir(pa):
        print('|    ' * de + '+--' + i)
        di = pa + '/' + i
        if isdir(di):
            fileTree(di, de + 1)


def filetree(pa=r"D:\User", de=0, spa=r"D:\Programming\User"):
    def fi(pa, de):
        from os.path import isdir
        from os import listdir
        if de == 0:
            f.write('文件夹:' + pa + '\n')
        for i in listdir(pa):
            f.write('|    ' * de + '+--' + i + '\n')
            di = pa + '/' + i
            if isdir(di):
                fi(di, de + 1)

    f = open(spa + '\\' + '_'.join(pa.replace(':',
             '').split(sep='\\')) + '_filetree2.txt', 'w')
    fi(pa, de)
    f.close()


def Gravitation(m1=1, m2=1, r=1, G=6.67e-11):
    return G * m1 * m2 / (r ** 2)


def ElectricForce(e1=1, e2=1, r=1, k=8.99e9):
    return k * e1 * e2 / (r ** 2)


def SurroundSpeed(m=1, r=1, f=1):
    return (f * r / m) ** 0.5


def ElectricFieldStrength(q=1, r=1, k=8.99e9):
    return q * k / (r * r)


def ElectronSpeed(u=1, e=1.6e-19, m=9.1e-31):
    return (2 * u * e / m) ** 0.5


if __name__ == '__main__':
    # print(gethtmlbypyppteer())
    print(math.log(2.7))
    print('{:e}'.format(ElectronSpeed(u=300)))

if __name__ == '__main__':
    hehe()