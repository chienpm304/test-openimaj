package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.line.Line2d;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.zip.DeflaterOutputStream;

public class LineHolder implements Comparable<LineHolder>{
    public List<Line2d> lines;
    public double rank = 0;
    public LineHolder(Line2d line1, Line2d line2, Line2d line3, Line2d line4){
        lines = new ArrayList<>();
        lines.add(line1);
        lines.add(line2);
        lines.add(line3);
        lines.add(line4);
    }

    @Override
    public int compareTo(LineHolder other) {
        return Double.compare(this.rank, other.rank);
    }

    @Override
    public String toString() {
        return "LineHolder{" +
                "lines=" + lines +
                ", rank=" + rank +
                '}';
    }
}
