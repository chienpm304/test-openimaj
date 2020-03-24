package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.processor.PixelProcessor;

public class ConvertAbsolute implements PixelProcessor<Float> {
    float alpha = 1.0f;
    float beta = 0f;

    public ConvertAbsolute(float alpha, float beta) {
        this.alpha = alpha;
        this.beta = beta;
    }

    @Override
    public Float processPixel(Float pixel) {
        float tmp = alpha * pixel + beta;
        tmp = tmp > 1f ? 1f : (tmp < 0 ? 0f : tmp);
        return tmp;
    }
}
