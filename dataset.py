import tensorflow as tf


def main():



    print(tf.__version__)
    dataset = tf.data.Dataset.from_tensor_slices([1, 2, 3])
    for element in dataset:
        print(element)
        tf.Tensor(1, shape=(), dtype=int32)
        tf.Tensor(2, shape=(), dtype=int32)
        tf.Tensor(3, shape=(), dtype=int32)






if __name__=="__main__":
    main()