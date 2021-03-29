# 一次元アーキテクチャにおける量子ビット割り当て問題

## 動作環境
- mac環境で動作確認を行っています
- windows環境やlinux環境での動作は保証しません
<br>

## 使っているライブラリ
- itertools
- queue
- amplify
- tkinter
<br>

## 実行方法
1. ファイル *application.py* の，`client.token`のコメントアウトを外してアクセストークンを入力する
<br>

2. *application.py* を実行する

```shell
$ python application.py
```
<br>

3. 読み込むサンプルを選択し，「読み込む」ボタンを押す

> *qasm/ex1.txt* : 2ビット加算器

> *qasm/ex2.txt* : grover探索による6ビット回文生成器
<br>

4. タイムアウト，制約の重みパラメータを調整し，「実行」ボタンを押す

> *qasm/ex2.txt* を選んだ場合，タイムアウトに加えて30秒~1分程度かかります．
<br>

5. 出力結果に満足がいかない場合，4に戻る
<br>

6. 満足のいく計算結果が得られたら，「保存」ボタンを押す

> *qasm/ex〇〇_output.txt* から確認できます
<br>

### ソルバーの切り替え
894~896行目のコメントアウトを切り替えることで，古典的解法や厳密解法も実行することができます

<br>

### レイヤーの表示

484行目の"False"を"True"に書き換えると，レイヤー間の仕切りも描画されます

<br>

## 実行結果

- 「読み込む」ボタンを押した結果

![image](https://user-images.githubusercontent.com/50867811/112040235-50a90780-8b88-11eb-807b-3cdd358b1b3a.png)

- パラメータ調整後，「実行」ボタンを押した結果

![image](https://user-images.githubusercontent.com/50867811/112040850-096f4680-8b89-11eb-9efa-a0e69831397c.png)

- 「保存」ボタンを押した結果

> SWAPゲートをCXゲートに展開した結果が描画されます

![image](https://user-images.githubusercontent.com/50867811/112040949-2572e800-8b89-11eb-97b8-6096ea563abe.png)

<br>

## 提出前チェック


- [x] README.mdの手順通りにして，プログラムが実行できる
- [x] 説明用スライドを用意した 
- [x] アクセストークンはリポジトリに含まれていない
- [x] MIT Licenseにした
