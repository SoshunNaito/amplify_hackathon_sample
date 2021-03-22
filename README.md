# 一次元アーキテクチャへの搭載に向けた量子回路設計支援アプリケーション

## 動作環境
- mac環境で動作確認を行っています
- windows環境やlinux環境での動作は保証しません

## 必要なライブラリ
- amplify
- tkinter

## 実行方法
1. ファイル *chapter5.py* の，`client.token`のコメントアウトを外してアクセストークンを入力する

2. *chapter5.py* を実行する

```shell
$ python chapter5.py
```

3. 読み込むサンプルを選択し，「読み込む」ボタンを押す

> *qasm/ex1.txt* : 2ビット加算器

> *qasm/ex2.txt* : grover探索による6ビット回文生成器

4. タイムアウト，制約の重みパラメータを調整し，「実行」ボタンを押す

> *qasm/ex2.txt* を選んだ場合，タイムアウトに加えて30秒~1分程度かかります．

5. 出力結果に満足がいかない場合，4に戻る

6. 満足のいく計算結果が得られたら，「保存」ボタンを押す

> *qasm/ex〇〇_output.txt* から確認できます

## 実行結果

- 「読み込む」ボタンを押した結果

![image](https://user-images.githubusercontent.com/50867811/112040235-50a90780-8b88-11eb-807b-3cdd358b1b3a.png)

- パラメータ調整後，「実行」ボタンを押した結果

![image](https://user-images.githubusercontent.com/50867811/112040850-096f4680-8b89-11eb-9efa-a0e69831397c.png)

- 「保存」ボタンを押した結果

> SWAPゲートをCXゲートに展開した結果が描画されます

![image](https://user-images.githubusercontent.com/50867811/112040949-2572e800-8b89-11eb-97b8-6096ea563abe.png)



## 提出前チェック


- [x] README.mdの手順通りにして，プログラムが実行できる
- [x] 説明用スライドを用意した 
- [ ] アクセストークンはリポジトリに含まれていない
- [x] MIT Licenseにした
