# 一次元アーキテクチャへの搭載に向けた量子回路設計支援アプリケーション

## 実行方法

1. ファイル *chapter5.py* の，`client.token`のコメントアウトを外してアクセストークンを入力する

2. *chapter5.py* を実行する

```shell
$ python chapter5.py
```

3. 読み込むサンプルを選択し，「読み込む」ボタンを押す
- *qasm/ex1.txt* : 2ビット加算器
- *qasm/ex2.txt* : grover探索による8ビット回文生成器

4. タイムアウト，制約の重みパラメータを調整し，「実行」ボタンを押す

5. 満足のいく計算結果が得られたら，「保存」ボタンを押す
- *qasm/ex〇〇_output.txt* に書き出される

## 実行結果

```
initial
2   5 | 1 3   |     4 
      |   4 8 |       
      |     7 |   2   
---------------------
  3 8 | 5     |   9 2 
      |   9   | 7     
      |       | 4 5   
---------------------
8 6   | 9 7   |       
9 5   |       |   3 1 
    4 |       |       
solved
2 8 5 | 1 3 9 | 6 7 4 
6 7 3 | 2 4 8 | 5 1 9 
4 1 9 | 6 5 7 | 3 2 8 
---------------------
7 3 8 | 5 6 4 | 1 9 2 
5 4 2 | 3 9 1 | 7 8 6 
1 9 6 | 7 8 2 | 4 5 3 
---------------------
8 6 1 | 9 7 3 | 2 4 5 
9 5 7 | 4 2 6 | 8 3 1 
3 2 4 | 8 1 5 | 9 6 7 
```

## 提出前チェック


- [ ] README.mdの手順通りにして，プログラムが実行できる
- [ ] 説明用スライドを用意した 
- [ ] アクセストークンはリポジトリに含まれていない
- [ ] MIT Licenseにした
