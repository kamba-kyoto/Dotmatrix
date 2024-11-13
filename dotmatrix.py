#椛 希代香　23261037
#モジュールのインポート
from Bio import SeqIO
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from io import StringIO # アップロードファイル操作用

#ウィンドウサイズの設定
#win = 10

#配列の読み込み
def dotmatrix(f1,f2,win):
    record1 = next(SeqIO.parse(f1, "fasta"))
    record2 = next(SeqIO.parse(f2, "fasta"))
    width = 500
    height = 500
    seq1 = record1.seq
    seq2 = record2.seq
    len1 = len(seq1) - win + 1
    len2 = len(seq2) - win + 1

    #ゼロ行列の作成
    image = np.zeros((height, width))

    #ハッシュの作成
    hash = {}
    for x in range(len1):
        subseq1 = str(seq1[x:x + win])
        if subseq1 not in hash:
            hash[subseq1] = []
        hash[subseq1].append(x)

    #部分配列の比較
    for y in range(len2):
        subseq2 = seq2[y:y + win]
        py = int(y/len2*height)
        if subseq2 in hash:
            for x in hash[subseq2]:
                px = int(x/len1*width)
                image[py, px] = 1

    #結果の表示
    plt.imshow(image, extent=(1,len1,len2,1), cmap="Grays")
    #plt.show() 最後の行をコメントアウト
    st.pyplot(plt) 
    #書き加える：Streamlit上にMatplotlibを表示

st.title("Dot matrix") # タイトル
# 配列ファイルのアップローダ
file1=st.sidebar.file_uploader("Sequence file 1:")
file2=st.sidebar.file_uploader("Sequence file 2:")
win=st.sidebar.slider("Window size:",4,100,10) # スライダ ー

if file1 and file2: #2つのファイルがアップロードされていれば
    with StringIO(file1.getvalue().decode("utf-8")) as f1,\
         StringIO(file2.getvalue().decode("utf-8")) as f2:
        dotmatrix(f1,f2,win) # 関数呼び出し