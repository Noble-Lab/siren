module load java/8u151

PARAM=$1
mzXML=$2

java -jar -Xmx8G ../java/DIA_Umpire_SE2.java $mzXML $PARAM 


