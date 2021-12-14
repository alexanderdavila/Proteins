import tensorflow as tf


def main():


    print("test")
    print(tf.__version__)
    print("test")
    dataset = tf.data.Dataset.from_tensor_slices([1, 2, 3])
    print(list(dataset.as_numpy_iterator()))






if __name__=="__main__":
    main()