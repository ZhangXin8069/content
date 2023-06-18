# -*- coding: utf-8 -*-

import os
from pdfminer.pdfparser import PDFParser, PDFDocument
from pdfminer.pdfinterp import PDFTextExtractionNotAllowed
from pdfminer.layout import *
from pdfminer.converter import PDFPageAggregator
from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
import time
import random
import hashlib
import sys
import importlib
importlib.reload(sys)


# **********翻译部分********************
def fanyi(query):
    import http.client
    import hashlib
    import urllib
    import random
    import json
    appid = ''  # !!!!补充
    secretKey = ''  # !!!!补充
    httpClient = None
    myurl = '/api/trans/vip/translate'
    q = query
    fromLang = 'auto'
    toLang = 'zh'
    salt = random.randint(32768, 65536)

    sign = appid + q + str(salt) + secretKey
    m1 = hashlib.md5()
    m1.update(sign.encode())
    sign = m1.hexdigest()
    myurl = myurl + '?appid=' + appid + '&q=' + urllib.parse.quote(
        q) + '&from=' + fromLang + '&to=' + toLang + '&salt=' + str(salt) + '&sign=' + sign
    # print(urllib.parse.quote(q))
    try:
        httpClient = http.client.HTTPConnection('api.fanyi.baidu.com')
        httpClient.request('GET', myurl)
        # response是HTTPResponse对象
        response = httpClient.getresponse()
        html = response.read()  # bytes
        # print("html:  ",type(html),html)
        html_str = html.decode()  # bytes to str
        # print("html_str:  ",type(html_str),html_str)
        html_dict = json.loads(html_str)  # str to dict
        # print("html_dist:  ",type(html_dict),html_str)
        # result_ori = html_dict["trans_result"][0]["src"]
        # result_tar = html_dict["trans_result"][0]["dst"]
        # print(html_dict["trans_result"])
        result_tar = ''
        for i in html_dict["trans_result"]:
            result_tar += i["dst"]
        # print(result_ori, " --> ", result_tar)
        print("翻译文本: " + result_tar)
        print("*" * 100)
        return result_tar
    except Exception as e:
        print(e)
        return ''
    finally:
        if httpClient:
            httpClient.close()


'''
解析pdf文件，获取文件中包含的各种对象
'''

# 解析pdf文件函数


def parse(pdf_path):
    textName = pdf_path.split('\\')[-1].split('.')[0] + '.txt'
    fp = open(pdf_path, 'rb')  # 以二进制读模式打开
    # 用文件对象来创建一个pdf文档分析器
    parser = PDFParser(fp)
    # 创建一个PDF文档
    doc = PDFDocument()
    # 连接分析器 与文档对象
    parser.set_document(doc)
    doc.set_parser(parser)

    # 提供初始化密码
    # 如果没有密码 就创建一个空的字符串
    doc.initialize()

    # 检测文档是否提供txt转换，不提供就忽略
    if not doc.is_extractable:
        raise PDFTextExtractionNotAllowed
    else:
        # 创建PDf 资源管理器 来管理共享资源
        rsrcmgr = PDFResourceManager()
        # 创建一个PDF设备对象
        laparams = LAParams()
        device = PDFPageAggregator(rsrcmgr, laparams=laparams)
        # 创建一个PDF解释器对象
        interpreter = PDFPageInterpreter(rsrcmgr, device)

        # 用来计数页面，图片，曲线，figure，水平文本框等对象的数量
        num_page, num_image, num_curve, num_figure, num_TextBoxHorizontal = 0, 0, 0, 0, 0

        # 循环遍历列表，每次处理一个page的内容
        for page in doc.get_pages():  # doc.get_pages() 获取page列表
            num_page += 1  # 页面增一
            print("\r\n>> 当前页：", num_page)
            interpreter.process_page(page)
            # 接受该页面的LTPage对象
            layout = device.get_result()
            for x in layout:
                if isinstance(x, LTImage):  # 图片对象
                    num_image += 1
                if isinstance(x, LTCurve):  # 曲线对象
                    num_curve += 1
                if isinstance(x, LTFigure):  # figure对象
                    num_figure += 1
                if isinstance(x, LTTextBoxHorizontal):  # 获取文本内容
                    num_TextBoxHorizontal += 1  # 水平文本框对象增一
                    results = x.get_text()
                    print(results.replace('\n', ''))
                    # 保存文本内容
                    with open(textName, 'a+', encoding='utf8') as f:
                        results = x.get_text()
                        f.write(results.replace('\n', '') + '\n')
        print('对象数量：\n', '页面数：%s\n' % num_page, '图片数：%s\n' % num_image, '曲线数：%s\n' % num_curve, '水平文本框：%s\n'
              % num_TextBoxHorizontal)


if __name__ == '__main__':

    pdf_path = r'sympy-docs-pdf-1.10.1.pdf'
    rootPath = '\\'.join(pdf_path.split('\\')[:-1]) if "\\" in pdf_path else ''
    textName = pdf_path.split('\\')[-1].split('.')[0] + '.txt'
    print(">> 当前文件：", os.path.join(rootPath, textName))
    if os.path.exists(os.path.join(rootPath, textName)):
        print(">> 删除：", textName)
        os.remove(os.path.join(rootPath, textName))
    if os.path.exists(os.path.join(rootPath, "translate.txt")):
        print(">> 删除：", "translate.txt")
        os.remove(os.path.join(rootPath, "translate.txt"))

    parse(pdf_path)

    with open(textName, 'r', encoding='utf8') as f:
        content = f.read()
        results = content.split('.')
        for i in results:
            res = fanyi(i)
            with open("translate.txt", 'a+', encoding='utf8') as fp:
                fp.write(res + '\n')
            time.sleep(1)
