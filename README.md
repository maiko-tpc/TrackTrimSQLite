#　TrackTrimSQLite
## 概要
Garfield++の飛跡生成クラスTrackSrimがうまく動作していないようだったので、
代わりになるTrackクラス (TrackTrimSQLite) をつくりました。
TrackTrimSQLiteではSQLiteデータベースのレコードをもとに
TrackSrimと同じような形式で飛跡を生成します。
飛跡生成の前に、TRIMの計算結果をSQLiteのデータベースに保存しておく必要があります。

## 説明
### makedb.cpp
TRIMを使って計算した、COLLISON.txtとRANGE_3D.txtの２つのファイルから
イオンと原子の衝突を記録するデータベースを作成するためのプログラムです。
実行時に、入力ファイルのディレクトリと出力するデータベースのファイル名を引数として指定します。
動作にはSQLiteのC言語APIが必要です。

### testTrackTrimSQLite.cpp
TrackTrimSQLiteの使用例です。
計算結果はROOTのTTree形式で出力します。
動作にはSQLiteのC言語のAPIの他に、ROOTとGarfield++が必要です。