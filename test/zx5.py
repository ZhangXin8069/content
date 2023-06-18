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


if __name__ == '__main__':
    pah = '/home/zhangxin/桌面'
    # absolute path
    lit = findAllfiles(pah=pah)
    fileSort(file_lit=lit, pah=pah)
    fileStatisticsbyPython_zx(pah=pah)
