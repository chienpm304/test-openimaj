package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;


public class Tetragram {//implements Comparable<Tetragram> {
    private Point2d tl = new Point2dImpl(),
            tr = new Point2dImpl(),
            br = new Point2dImpl(),
            bl = new Point2dImpl();

    public Tetragram(Point2d topLeft, Point2d topRight,
                     Point2d bottomRight, Point2d bottomLeft) {
        this.tl = topLeft;
        this.tr = topRight;
        this.br = bottomRight;
        this.bl = bottomLeft;

    }

    public Tetragram(List<Point2d> points){
        if(points.size() !=4)
            throw new RuntimeException("Points passed must has 4 elements");
        tl = points.get(0);
        tr = points.get(1);
        br = points.get(2);
        bl = points.get(3);
    }

    public Tetragram(float[][] data){
        try{
            tl.setX(data[0][0]);
            tl.setY(data[0][1]);
            tr.setX(data[1][0]);
            tr.setY(data[1][1]);
            br.setX(data[2][0]);
            br.setY(data[2][1]);
            bl.setX(data[3][0]);
            bl.setY(data[3][1]);

        }catch (Exception e){
            e.printStackTrace();
            throw new RuntimeException("Data array must be shape of float[4][2]");
        }
    }

    public List<Point2d> toList(){
        List<Point2d> res = new ArrayList<>();
        res.add(tl);
        res.add(tr);
        res.add(br);
        res.add(bl);
        return res;
    }

    public float[][] toArray(){
        float [][]arr = new float[4][2];
        arr[0][0] = tl.getX();
        arr[0][1] = tl.getY();

        arr[1][0] = tr.getX();
        arr[1][1] = tr.getY();

        arr[2][0] = br.getX();
        arr[2][1] = br.getY();

        arr[3][0] = bl.getX();
        arr[3][1] = bl.getY();
        return arr;
    }

    @Override
    public String toString() {
        return "Tetragram{" +
                "tl=" + tl +
                ", tr=" + tr +
                ", br=" + br +
                ", bl=" + bl +
                '}';
    }

    public Point2d getTopLeft() {
        return tl;
    }

    public Point2d getTopRight() {
        return tr;
    }

    public Point2d getBottomRight() {
        return br;
    }

    public Point2d getBottomLeft() {
        return bl;
    }

    public Tetragram scale(float scaleFactor) {
        float[][] tmp = toArray();
        for(int i =0; i < 4; i++){
            for(int j = 0; j < 2; j++)
                tmp[i][j] *= scaleFactor;
        }
        return new Tetragram(tmp);
    }

    /**
     * @return a collection of lines in order: top, right, bottom, left
     */
    public List<Line2d> toLineList() {
        ArrayList<Line2d> lines = new ArrayList<Line2d>();
        lines.add(new Line2d(tl, tr));
        lines.add(new Line2d(tr, br));
        lines.add(new Line2d(br, bl));
        lines.add(new Line2d(bl, tl));
        return lines;
    }

//    @Override
//    public int compareTo(Tetragram o) {
//        Polygon thizPoly = new Polygon(this.toList());
//        Polygon thatPoly = new Polygon(o.toList());
//        return Double.compare(thizPoly.calculateArea(), thatPoly.calculateArea());
//    }

    public static Tetragram createRectangleBounding(int width, int height) {
        float[][] rect = new float[][]{
                {0, 0},
                {width-1, 0},
                {width-1, height-1},
                {0, height-1}
        };

        return new Tetragram(rect);
    }
}
