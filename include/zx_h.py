# import math
import nest_asyncio
nest_asyncio.apply()


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

# def filetree(pa=r"./", de=0, spa=r"./"):
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

#     f = open(spa + '\\' + '_'.join(pa.replace(':', '').split(sep='\\')) + '_filetree2.xlsx', 'w')
#     fi(pa, de)
#     f.close()
# # filetree()
# def fileTree(pa=r"./", de=0):
#     from os.path import isdir
#     from os import listdir
#     if de == 0:
#         print('文件夹:' + pa)
#     for i in listdir(pa):
#         print('|    ' * de + '+--' + i)
#         di = pa + '/' + i
#         if isdir(di):
#             fileTree(di, de + 1)
# print(filetree())


def gethtmlbyrequests(
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
#     from pyautogui import alert
    # pa = r"C:/Users/{}/Desktop/{}".format(getuser(), wo)
    # pa = r"D:/Desktop/{}".format(wo)
    pa = "./{}".format(wo)
    try:
        mkdir(pa)
    except FileExistsError:
        pass
    # li = findall(r'\"url\":\"(https://[^"]*?\.[pngje]{3,4}[^"]*?)\"', html)
    li = findall(r'\"(https://[^"]*\.[pngje]{3,4}[^"]*)\"', html)
#     alert('已找到{}张图片'.format(len(li)))
    print('已找到{}张图片'.format(len(li)))

    for i in range(nu):
        try:
            pb = get(li[i]).content
            path = pa+r'/{}_{}_{}.png'.format(wo, i + 1, nu)
            f = open(path, 'wb')
            f.write(pb)
            f.close()
            print('downloading...'+str(path)+'\nurl:'+str(li[i]))
        except Exception as e:
            print(e)
            continue
#     alert('下载完毕(可能漏下)')
    print('下载完毕(可能漏下)')


def gethtmlbypyppteer(
        url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    import asyncio
    import nest_asyncio
    nest_asyncio.apply()
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
    pa = r'D:\User'
    # pa = r'D:\DeskTop\2021级核物二班 活动记录'
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


def FileStatisticsbyPython_zx():
    from os import listdir, path, chdir

    from pandas import DataFrame

    # pa = str(input('postion(\\):'))
    pa = r'D:\Monitor'
    pa = r'D:\DeskTop\2021级核物二班 活动记录'
    # pa = r"D:\Monitor - 副本"
    pa0 = pa.split(sep=('\\'))[-1]
    # print(pa0)
    chdir(pa)
    li = []
    for i in listdir(pa):
        if path.isdir(pa + '\\' + i):
            li.append(i)
    lc = []
    for i in li:
        lt = list('=HYPERLINK(\".\\'+pa0+'\\' + i + '\\' + j + '\",\"' + j + '\")' for j in
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
    df.to_excel('..\\' + 'fileStatisticsbyPython_zx.xlsx', na_rep=' ')
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


def fileTree(pa=r"./", de=0):
    from os.path import isdir
    from os import listdir
    if de == 0:
        print('文件夹:' + pa)
    for i in listdir(pa):
        print('|    ' * de + '+--' + i)
        di = pa + '/' + i
        if isdir(di):
            fileTree(di, de + 1)


def filetree(pa=r"./", de=0, spa=r"./"):
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


def Unm(a=1, b=10, n=10000):
    from random import uniform
    x = []
    for i in range(0, n):
        x.append(uniform(a, b))
    return x


def pufg(a=1, l=1, n=1000):
    from random import uniform
    from math import pi, cos
    t = 0
    for i in range(0, n):
        if l*cos(uniform(0, pi/2)) >= uniform(0, a):
            t += 1
    return 2*n/t


def cir_s(n=10000):
    from random import uniform
    from math import pi, cos
    t = 0
    for i in range(0, n):
        if(uniform(-1, 1)**2+uniform(-1, 1)**2) <= 1:
            t += 1
    return t/n


def fileStatisticsbyPython_zx(pah='./'):
    from os import listdir
    from os.path import isdir, join, split
    from pandas import DataFrame
    lit = []
    for i in listdir(pah):
        if isdir(join(pah, i)):
            lit.append(join(pah, i))
    print(lit)
    lc = []
    for i in lit:
        lt = list('=HYPERLINK(\"' + join(pah, i, j) +
                  '\",\"' + j + '\")' for j in listdir(join(pah, i)))
        lt.insert(0, str(len(listdir(join(pah, i)))))
        lc.append(lt)
    df = DataFrame(data=lc, index=lit)
    df.to_excel(join(pah, 'fileStatisticsbyPython_zx.xlsx'), na_rep=' ')
    print(df)


def findAllfiles(pah):
    from os import walk, path
    lit = []
    for foldName, subfolders, filenames in walk(pah):
        for filename in filenames:
            try:
                lit.append(path.join(foldName, filename))
            except Exception as e:
                print(e)
                pass
    return lit


def fileSort(pah=[], want_lit=['pdf', 'pptx'], file_lit=[]):
    from os import mkdir, rename
    from os.path import join, split
    for i in want_lit:
        try:
            i = mkdir(join(pah, i))
        except Exception as e:
            print(e)
        finally:
            for j in file_lit:
                (_, j_) = split(j)
                if i in j_:
                    rename(j, join(pah, i, j_))


# if __name__ == '__main__':
#     pah = '/home/zhangxin/桌面'
#     # absolute path
#     lit = findAllfiles(pah=pah)
#     fileSort(file_lit=lit, pah=pah)
#     fileStatisticsbyPython_zx(pah=pah)

if __name__ == '__main__':
    # print(Unm(), len(Unm()))
    # print(pufg(n=10000))
    # fileStatisticsbyPython_zx()
    filetree(pa='D:\\user')
    # FileStatisticsbyPython_zx()
    # filetree()
# =============================================================================
#     print(gethtmlbypyppteer())
#     print(math.log(2.7))
#     # print('{:e}'.format(ElectronSpeed(u=300)))
# =============================================================================
    pass
