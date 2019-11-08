sudo apt-get install imagemagick

rm -rf car3div10
rm -rf car3div4
rm -rf car3div2
rm -rf car3x2
rm -rf car3x3
rm -rf car3x4
rm -rf car3x5
rm -rf car3x10

mkdir car3div10
mkdir car3div4
mkdir car3div2
mkdir car3x2
mkdir car3x3
mkdir car3x4
mkdir car3x5
mkdir car3x10


cp -r car3/*.pgm  car3div10
cp -r car3/*.pgm  car3div4
cp -r car3/*.pgm  car3div2
cp -r car3/*.pgm  car3x2
cp -r car3/*.pgm  car3x3
cp -r car3/*.pgm  car3x4
cp -r car3/*.pgm  car3x5
cp -r car3/*.pgm  car3x10

sudo mogrify -resize 10% car3div10/*.pgm
sudo mogrify -resize 25% car3div4/*.pgm
sudo mogrify -resize 50% car3div2/*.pgm
sudo mogrify -resize 200% car3x2/*.pgm
sudo mogrify -resize 300% car3x3/*.pgm
sudo mogrify -resize 400% car3x4/*.pgm
sudo mogrify -resize 500% car3x5/*.pgm
sudo mogrify -resize 1000% car3x10/*.pgm