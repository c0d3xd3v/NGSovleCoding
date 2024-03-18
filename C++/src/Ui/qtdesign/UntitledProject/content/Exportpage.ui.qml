import QtQuick 6.6
import QtQuick.Controls 6.6
import QtQuick.Layouts 6.6

Item {
    id: exportPage
    visible: true
    BackgroundToolbar {
        id: exportToolbar
        height: 48
        visible: true
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.leftMargin: 0
        anchors.rightMargin: 0
        anchors.topMargin: 0
        contentHeight: 37

        RowLayout {
            id: exportRowLayout
            visible: true
            anchors.fill: parent

            ToolButton {
                id: exportBackButton
                width: 42
                height: 42
                visible: true
                text: qsTr("Tool Button")
                icon.color: "#000000"
                Layout.fillHeight: false
                Layout.fillWidth: false
                padding: 1
                spacing: 6
                icon.height: 36
                icon.width: 36
                icon.source: "/home/kai/Development/github/NGSovleCoding/C++/src/Ui/icons/arrow_back.svg"
                display: AbstractButton.IconOnly
                onClicked: swipeView.currentIndex = 0
            }

            Item {
                width: 200
                height: 200
                Layout.fillHeight: true
                Layout.fillWidth: true
            }
        }
    }

    Pane {
        id: exportMainpage
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: exportToolbar.bottom
        anchors.bottom: parent.bottom
        anchors.topMargin: 0
    }
}
