#!/bin/bash


infile=$1

skip_begin=571
skip_end=2
tail_start=$((skip_begin+1))

tail --lines=+$tail_start $infile | head --lines=-$skip_end | awk '

	function max(x,y){

		if(x>y) return x
		else	return y
	}

	{
		totalHb  += $1  + 2.0*( $2+$14+$17+$18+$20+$24+$27+$28+$30+$34+$37+$38+$40)
		totalKr  += $11 + 2.0*( $4+ $7+ $8+$10+$12+$25+$27+$29+$30+$35+$37+$39+$40)
		totalKni += $21 + 2.0*( $5+ $7+ $9+$10+$15+$17+$19+$20+$22+$36+$38+$39+$40)
		totalGt  += $31 + 2.0*( $6+ $8+ $9+$10+$16+$18+$19+$20+$26+$28+$29+$30+$32)
	}
	END{

		totalSum = totalHb + totalKr + totalKni + totalGt
		lambda = (totalHb + max(totalKr, totalGt)) / totalSum
		print lambda, totalHb, totalKr, totalKni, totalGt

	}
'

